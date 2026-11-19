//! Random access to local UCSC twoBit v0 references.
//!
//! Coordinates are zero-based and half-open, and contig names are literal.
//! The file carries no assembly identifier. Verify the file's digest against a
//! trusted reference receipt before treating its coordinates as assembly-specific.
use anyhow::{ensure, Context, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, Read, Seek, SeekFrom},
    ops::Range,
    path::Path,
};

const SIGNATURE: u32 = 0x1A412743;
const BUFFER_BYTES: usize = 64 * 1024;
const BASES: [u8; 4] = *b"TCAG";

const fn base_index(byte: u8, offset: usize) -> usize {
    ((byte >> (6 - 2 * offset)) & 3) as usize
}

const fn decode_table() -> [[u8; 4]; 256] {
    let mut table = [[0; 4]; 256];
    let mut byte = 0;
    while byte < table.len() {
        let mut position = 0;
        while position < 4 {
            table[byte][position] = BASES[base_index(byte as u8, position)];
            position += 1;
        }
        byte += 1;
    }
    table
}

const fn count_table() -> [[u8; 4]; 256] {
    let mut table = [[0; 4]; 256];
    let mut byte = 0;
    while byte < table.len() {
        let mut position = 0;
        while position < 4 {
            table[byte][base_index(byte as u8, position)] += 1;
            position += 1;
        }
        byte += 1;
    }
    table
}

const DECODE: [[u8; 4]; 256] = decode_table();
const COUNTS: [[u8; 4]; 256] = count_table();

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum ByteOrder {
    Big,
    Little,
}

impl ByteOrder {
    fn read_u32(self, source: &mut impl Read) -> Result<u32> {
        let mut bytes = [0; 4];
        source.read_exact(&mut bytes)?;
        Ok(match self {
            Self::Big => u32::from_be_bytes(bytes),
            Self::Little => u32::from_le_bytes(bytes),
        })
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Contig<'a> {
    pub name: &'a str,
    pub length: u64,
    pub n_blocks: usize,
    pub mask_blocks: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct BaseCounts {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
}

struct Block {
    start: u32,
    end: u32,
}

struct Record {
    length: u32,
    ns: Vec<Block>,
    masks: Vec<Block>,
    packed_offset: u64,
}

/// Validated index and mask tables over a seekable reference.
/// Packed DNA is read on demand using a reusable bounded buffer.
pub struct TwoBitReference<R> {
    source: R,
    file_bytes: u64,
    order: ByteOrder,
    records: BTreeMap<String, Record>,
    scratch: Vec<u8>,
}

impl TwoBitReference<File> {
    pub fn open(path: &Path) -> Result<Self> {
        let source = File::open(path).with_context(|| format!("open {}", path.display()))?;
        Self::new(source).with_context(|| format!("read twoBit file {}", path.display()))
    }
}

impl<R: Read + Seek> TwoBitReference<R> {
    /// Validate the file header, index and block tables, without unpacking DNA.
    pub fn new(mut source: R) -> Result<Self> {
        let file_bytes = source.seek(SeekFrom::End(0))?;
        let (order, records) = parse(&mut BufReader::new(&mut source), file_bytes)?;
        Ok(Self {
            source,
            file_bytes,
            order,
            records,
            scratch: vec![0; BUFFER_BYTES],
        })
    }

    pub fn file_bytes(&self) -> u64 {
        self.file_bytes
    }

    pub fn byte_order(&self) -> ByteOrder {
        self.order
    }

    /// Enumerate sequences in ascending literal name order.
    pub fn contigs(&self) -> impl Iterator<Item = Contig<'_>> {
        self.records.iter().map(|(name, record)| Contig {
            name,
            length: u64::from(record.length),
            n_blocks: record.ns.len(),
            mask_blocks: record.masks.len(),
        })
    }

    /// Hash the complete source bytes. Call once when binding a reference receipt.
    pub fn sha256(&mut self) -> Result<String> {
        self.source.seek(SeekFrom::Start(0))?;
        let mut remaining = self.file_bytes;
        let mut digest = Sha256::new();
        while remaining != 0 {
            let bytes = remaining.min(self.scratch.len() as u64) as usize;
            let chunk = &mut self.scratch[..bytes];
            self.source.read_exact(chunk)?;
            digest.update(chunk);
            remaining -= bytes as u64;
        }
        ensure!(
            self.source.seek(SeekFrom::End(0))? == self.file_bytes,
            "twoBit file length changed while hashing"
        );
        Ok(format!("{:x}", digest.finalize()))
    }

    /// Decode an interval, applying N blocks before optional lowercase masking.
    pub fn fetch(
        &mut self,
        name: &str,
        start: u64,
        end: u64,
        strand: Strand,
        masked: bool,
    ) -> Result<String> {
        let record = find_record(&self.records, name)?;
        ensure!(
            start <= end && end <= u64::from(record.length),
            "twoBit interval must satisfy start <= end <= contig length"
        );
        let length = usize::try_from(end - start)?;
        if length == 0 {
            return Ok(String::new());
        }

        let mut output = Vec::new();
        output.try_reserve_exact(length)?;
        read_packed(
            &mut self.source,
            &mut self.scratch,
            record.packed_offset,
            start..end,
            |chunk, _, span| append_bases(chunk, span, &mut output),
        )?;
        for block in intersections(&record.ns, start, end) {
            output[(block.start - start) as usize..(block.end - start) as usize].fill(b'N');
        }
        if masked {
            for block in intersections(&record.masks, start, end) {
                output[(block.start - start) as usize..(block.end - start) as usize]
                    .make_ascii_lowercase();
            }
        }
        if strand == Strand::Reverse {
            output.reverse();
            for base in &mut output {
                *base = match *base {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    b'a' => b't',
                    b'c' => b'g',
                    b'g' => b'c',
                    b't' => b'a',
                    b'N' | b'n' => *base,
                    _ => unreachable!("twoBit decoder emits ACGTN only"),
                };
            }
        }
        Ok(String::from_utf8(output).expect("twoBit decoder emits ASCII"))
    }

    /// Count bases directly from packed bytes, with N blocks overriding the DNA.
    /// Masking does not affect base counts.
    pub fn base_counts(&mut self, name: &str) -> Result<BaseCounts> {
        let record = find_record(&self.records, name)?;
        let mut known = [0; 4];
        let mut unknown = 0;
        read_packed(
            &mut self.source,
            &mut self.scratch,
            record.packed_offset,
            0..u64::from(record.length),
            |chunk, base, span| {
                let mut position = span.start;
                for block in intersections(&record.ns, base, base + span.end as u64) {
                    let start = (block.start - base) as usize;
                    tally_span(chunk, position..start, &mut known);
                    position = (block.end - base) as usize;
                    unknown += block.end - block.start;
                }
                tally_span(chunk, position..span.end, &mut known);
            },
        )?;
        let [t, c, a, g] = known;
        Ok(BaseCounts {
            a,
            c,
            g,
            t,
            n: unknown,
        })
    }
}

fn find_record<'a>(records: &'a BTreeMap<String, Record>, name: &str) -> Result<&'a Record> {
    records
        .get(name)
        .with_context(|| format!("contig is absent from the twoBit index: {name}"))
}

fn parse(
    source: &mut (impl Read + Seek),
    file_bytes: u64,
) -> Result<(ByteOrder, BTreeMap<String, Record>)> {
    source.seek(SeekFrom::Start(0))?;
    let (order, index, index_end) = read_index(&mut *source, file_bytes)?;
    let mut physical_order = index.into_iter().peekable();
    let mut records = BTreeMap::new();
    while let Some((offset, name)) = physical_order.next() {
        let limit = physical_order.peek().map_or(file_bytes, |next| next.0);
        ensure!(
            index_end <= offset && offset < limit && limit <= file_bytes,
            "twoBit record offsets are duplicated or outside the record region"
        );
        source.seek(SeekFrom::Start(offset))?;
        let record = read_record((&mut *source).take(limit - offset), order, limit)
            .with_context(|| format!("read twoBit record {name}"))?;
        ensure!(
            records.insert(name, record).is_none(),
            "duplicate twoBit contig name"
        );
    }
    Ok((order, records))
}

fn read_index(source: impl Read, file_bytes: u64) -> Result<(ByteOrder, Vec<(u64, String)>, u64)> {
    let mut index_source = source.take(file_bytes);
    let mut signature = [0; 4];
    index_source.read_exact(&mut signature)?;
    let order = if signature == SIGNATURE.to_be_bytes() {
        ByteOrder::Big
    } else {
        ensure!(
            signature == SIGNATURE.to_le_bytes(),
            "twoBit signature mismatch"
        );
        ByteOrder::Little
    };
    ensure!(
        order.read_u32(&mut index_source)? == 0,
        "only twoBit version 0 is supported"
    );
    let count = order.read_u32(&mut index_source)?;
    ensure!(
        order.read_u32(&mut index_source)? == 0,
        "twoBit header reserved field must be zero"
    );
    ensure!(
        u64::from(count) <= index_source.limit() / 6,
        "twoBit sequence count exceeds index capacity"
    );

    let mut index = Vec::new();
    let mut name = [0; u8::MAX as usize];
    for _ in 0..count {
        let mut size = [0];
        index_source.read_exact(&mut size)?;
        let name = &mut name[..usize::from(size[0])];
        index_source.read_exact(name)?;
        ensure!(
            !name.is_empty() && name.iter().all(u8::is_ascii_graphic),
            "twoBit contig names must be nonempty printable nonwhitespace ASCII"
        );
        let name = std::str::from_utf8(name)?.to_owned();
        let offset = u64::from(order.read_u32(&mut index_source)?);
        index.push((offset, name));
    }
    index.sort_unstable_by_key(|(offset, _)| *offset);
    Ok((order, index, file_bytes - index_source.limit()))
}

fn read_record(
    mut source: std::io::Take<impl Read>,
    order: ByteOrder,
    limit: u64,
) -> Result<Record> {
    let length = order.read_u32(&mut source)?;
    let ns = read_blocks(&mut source, order, length, "N")?;
    let masks = read_blocks(&mut source, order, length, "mask")?;
    ensure!(
        order.read_u32(&mut source)? == 0,
        "twoBit record reserved field must be zero"
    );
    let remaining = source.limit();
    ensure!(
        u64::from(length).div_ceil(4) <= remaining,
        "packed DNA exceeds its twoBit record boundary"
    );
    Ok(Record {
        length,
        ns,
        masks,
        packed_offset: limit - remaining,
    })
}

fn read_blocks(
    source: &mut std::io::Take<impl Read>,
    order: ByteOrder,
    length: u32,
    kind: &str,
) -> Result<Vec<Block>> {
    let count = order.read_u32(source)?;
    ensure!(
        count <= length && u64::from(count) * 8 <= source.limit(),
        "{kind} block table exceeds its twoBit record boundary"
    );
    let mut blocks = Vec::new();
    blocks.try_reserve_exact(usize::try_from(count)?)?;
    for _ in 0..count {
        blocks.push(Block {
            start: order.read_u32(source)?,
            end: 0,
        });
    }
    let mut previous_end = 0;
    for block in &mut blocks {
        let size = order.read_u32(source)?;
        block.end = block
            .start
            .checked_add(size)
            .context("twoBit block overflows")?;
        ensure!(
            size != 0 && block.start >= previous_end && block.end <= length,
            "{kind} blocks must be positive, sorted, disjoint and in bounds"
        );
        previous_end = block.end;
    }
    Ok(blocks)
}

fn intersections(blocks: &[Block], start: u64, end: u64) -> impl Iterator<Item = Range<u64>> + '_ {
    let first = blocks.partition_point(|block| u64::from(block.end) <= start);
    blocks[first..]
        .iter()
        .take_while(move |block| u64::from(block.start) < end)
        .map(move |block| start.max(u64::from(block.start))..end.min(u64::from(block.end)))
}

/// Read only the packed bytes touched by an interval. Each callback receives
/// the coordinate of its first packed base and the eligible range in the chunk.
fn read_packed(
    source: &mut (impl Read + Seek),
    scratch: &mut [u8],
    packed_offset: u64,
    range: Range<u64>,
    mut consume: impl FnMut(&[u8], u64, Range<usize>),
) -> Result<()> {
    if range.is_empty() {
        return Ok(());
    }
    let mut byte_position = range.start / 4;
    let end_byte = range.end.div_ceil(4);
    source.seek(SeekFrom::Start(packed_offset + byte_position))?;
    while byte_position < end_byte {
        let bytes = (end_byte - byte_position).min(scratch.len() as u64) as usize;
        let chunk = &mut scratch[..bytes];
        source.read_exact(chunk)?;
        let base = byte_position * 4;
        let start = range.start.saturating_sub(base) as usize;
        let end = (range.end - base).min(4 * bytes as u64) as usize;
        consume(chunk, base, start..end);
        byte_position += bytes as u64;
    }
    Ok(())
}

fn append_bases(mut packed: &[u8], span: Range<usize>, output: &mut Vec<u8>) {
    let mut length = span.len();
    if span.start != 0 {
        let first = (4 - span.start).min(length);
        output.extend_from_slice(&DECODE[packed[0] as usize][span.start..span.start + first]);
        packed = &packed[1..];
        length -= first;
    }
    let full = length / 4;
    for &byte in &packed[..full] {
        output.extend_from_slice(&DECODE[byte as usize]);
    }
    let tail = length % 4;
    if tail != 0 {
        output.extend_from_slice(&DECODE[packed[full] as usize][..tail]);
    }
}

fn tally_span(packed: &[u8], span: Range<usize>, totals: &mut [u64; 4]) {
    let mut base = span.start;
    while base < span.end && !base.is_multiple_of(4) {
        totals[base_index(packed[base / 4], base % 4)] += 1;
        base += 1;
    }
    let complete = span.end / 4;
    if base / 4 < complete {
        for &byte in &packed[base / 4..complete] {
            for (total, count) in totals.iter_mut().zip(COUNTS[byte as usize]) {
                *total += u64::from(count);
            }
        }
        base = complete * 4;
    }
    while base < span.end {
        totals[base_index(packed[base / 4], base % 4)] += 1;
        base += 1;
    }
}
