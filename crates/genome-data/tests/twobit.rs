use genome_data::twobit::{BaseCounts, ByteOrder, Strand, TwoBitReference};
use std::{
    cell::Cell,
    io::{self, Cursor, Read, Seek, SeekFrom},
    rc::Rc,
};

const B_OFFSET: usize = 34;
const Z_OFFSET: usize = 52;
const A_OFFSET: usize = 68;

#[derive(Clone, Copy)]
enum FixtureOrder {
    Big,
    Little,
}

impl FixtureOrder {
    fn word(self, value: u32) -> [u8; 4] {
        match self {
            Self::Big => value.to_be_bytes(),
            Self::Little => value.to_le_bytes(),
        }
    }
}

fn fixture(order: FixtureOrder) -> Vec<u8> {
    let mut raw = Vec::new();
    for value in [0x1A41_2743, 0, 3, 0] {
        raw.extend(order.word(value));
    }
    // Index, physical record, and alphabetical orders all differ.
    for (name, offset) in [(b'z', Z_OFFSET), (b'a', A_OFFSET), (b'b', B_OFFSET)] {
        raw.extend([1, name]);
        raw.extend(order.word(offset as u32));
    }
    assert_eq!(raw.len(), B_OFFSET);
    for value in [5, 0, 0, 0] {
        raw.extend(order.word(value));
    }
    // AAAAA with nonzero padding in the last packed byte.
    raw.extend([0xaa, 0xbf]);
    assert_eq!(raw.len(), Z_OFFSET);
    for value in [0, 0, 0, 0] {
        raw.extend(order.word(value));
    }
    assert_eq!(raw.len(), A_OFFSET);
    // TCAGACGTG with N blocks [2,5), [8,9) and mask blocks [1,5), [7,9).
    for value in [9, 2, 2, 8, 3, 1, 2, 1, 7, 4, 2, 0] {
        raw.extend(order.word(value));
    }
    raw.extend([0x1b, 0x9c, 0xc0]);
    raw
}

#[test]
fn every_small_interval_matches_both_endian_encodings() {
    let uppercase = "TCNNNCGTN";
    let masked = "TcnnnCGtn";
    for order in [FixtureOrder::Big, FixtureOrder::Little] {
        let raw = fixture(order);
        let mut reference = TwoBitReference::new(Cursor::new(&raw)).unwrap();
        assert_eq!(reference.file_bytes(), raw.len() as u64);
        assert!(matches!(
            (order, reference.byte_order()),
            (FixtureOrder::Big, ByteOrder::Big) | (FixtureOrder::Little, ByteOrder::Little)
        ));
        assert_eq!(
            reference
                .contigs()
                .map(|c| (c.name, c.length, c.n_blocks, c.mask_blocks))
                .collect::<Vec<_>>(),
            vec![("a", 9, 2, 2), ("b", 5, 0, 0), ("z", 0, 0, 0)]
        );
        for start in 0..=uppercase.len() {
            for end in start..=uppercase.len() {
                for (soft_mask, expected) in [(false, uppercase), (true, masked)] {
                    assert_eq!(
                        reference
                            .fetch("a", start as u64, end as u64, Strand::Forward, soft_mask)
                            .unwrap(),
                        &expected[start..end]
                    );
                }
            }
        }
        assert_eq!(
            reference.fetch("a", 0, 9, Strand::Reverse, true).unwrap(),
            "naCGnnngA"
        );
        let empty = reference.fetch("z", 0, 0, Strand::Reverse, false).unwrap();
        assert!(empty.is_empty());
    }
}

#[test]
fn counts_exclude_n_overlays_and_unused_packed_bits() {
    let raw = fixture(FixtureOrder::Little);
    let mut reference = TwoBitReference::new(Cursor::new(raw)).unwrap();
    let counts_by_contig = [
        ("a", [0, 2, 1, 2, 4]),
        ("b", [5, 0, 0, 0, 0]),
        ("z", [0; 5]),
    ];
    for (name, expected) in counts_by_contig {
        let counts = reference.base_counts(name).unwrap();
        assert_eq!([counts.a, counts.c, counts.g, counts.t, counts.n], expected);
    }
    assert!(reference.base_counts("A").is_err());
    let invalid_intervals = [
        ("chrA", 0, 1),
        ("a", 2, 1),
        ("a", 0, 10),
        ("a", 0, u64::MAX),
    ];
    for (name, start, end) in invalid_intervals {
        let result = reference.fetch(name, start, end, Strand::Forward, false);
        assert!(result.is_err());
    }
}

#[test]
fn rejects_malformed_tables_and_truncated_payload() {
    let raw = fixture(FixtureOrder::Big);
    for (offset, value) in [
        (0, 0),
        (4, 1),
        (8, u32::MAX),
        (12, 1),
        (18, 16),
        (18, (B_OFFSET + 17) as u32),
        (18, raw.len() as u32),
        (24, Z_OFFSET as u32),
        (A_OFFSET + 4, u32::MAX),
        (A_OFFSET + 12, 1),
        (A_OFFSET + 16, 0),
        (A_OFFSET + 36, 7),
        (A_OFFSET + 44, 1),
    ] {
        let mut corrupt = raw.clone();
        corrupt[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
        assert!(
            TwoBitReference::new(Cursor::new(corrupt)).is_err(),
            "unexpected acceptance of {value} at offset {offset}"
        );
    }
    let mut duplicate_name = raw.clone();
    duplicate_name[23] = b'z';
    assert!(TwoBitReference::new(Cursor::new(duplicate_name)).is_err());
    assert!(TwoBitReference::new(Cursor::new(&raw[..raw.len() - 1])).is_err());
}

struct TracedReader {
    inner: Cursor<Vec<u8>>,
    bytes_read: Rc<Cell<usize>>,
    seeks: Rc<Cell<usize>>,
}

impl Read for TracedReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let count = self.inner.read(buf)?;
        self.bytes_read.set(self.bytes_read.get() + count);
        Ok(count)
    }
}

impl Seek for TracedReader {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.seeks.set(self.seeks.get() + 1);
        self.inner.seek(position)
    }
}

#[test]
fn short_queries_read_only_their_packed_span() {
    let bytes_read = Rc::new(Cell::new(0));
    let mut reference = TwoBitReference::new(TracedReader {
        inner: Cursor::new(fixture(FixtureOrder::Big)),
        bytes_read: Rc::clone(&bytes_read),
        seeks: Rc::default(),
    })
    .unwrap();
    bytes_read.set(0);
    let sequence = reference.fetch("a", 3, 6, Strand::Forward, false).unwrap();
    assert_eq!(sequence, "NNC");
    assert!(bytes_read.get() <= 2);
}

#[test]
fn partitioned_counts_match_scalar_oracle_across_packed_chunks() {
    let edge = 4 * 64 * 1024u64;
    let length = 2 * edge + 5;
    let mut raw = [0x1a41_2743u32, 0, 1, 0]
        .into_iter()
        .flat_map(u32::to_le_bytes)
        .collect::<Vec<_>>();
    raw.extend([1, b'x']);
    let n_blocks = [edge - 5..edge + 9, 2 * edge - 2..2 * edge + 2];
    raw.extend(22u32.to_le_bytes());
    raw.extend((length as u32).to_le_bytes());
    raw.extend(2u32.to_le_bytes());
    for block in &n_blocks {
        raw.extend((block.start as u32).to_le_bytes());
    }
    for block in &n_blocks {
        raw.extend(((block.end - block.start) as u32).to_le_bytes());
    }
    // The soft mask intersects an N block. Packed 0x1b spells TCAG and the
    // final byte has three nonzero padding symbols.
    for word in [1, edge as u32 - 4, 16, 0] {
        raw.extend(word.to_le_bytes());
    }
    let packed_bytes = length.div_ceil(4) as usize;
    raw.resize(raw.len() + packed_bytes, 0x1b);

    let runs = [
        (3, 0),
        (10, 0),
        (edge - 3, 1),
        (edge + 4, 2),
        (edge + 5, 2),
        (2 * edge - 1, 0),
        (2 * edge + 4, 2),
        (length, 1),
    ];
    let mut expected = [BaseCounts::default(); 4];
    let mut start = 0;
    for &(end, class) in &runs {
        let counts = &mut expected[class];
        for position in start..end {
            if n_blocks.iter().any(|block| block.contains(&position)) {
                counts.n += 1;
            } else {
                match position % 4 {
                    0 => counts.t += 1,
                    1 => counts.c += 1,
                    2 => counts.a += 1,
                    _ => counts.g += 1,
                }
            }
        }
        start = end;
    }

    let bytes_read = Rc::new(Cell::new(0));
    let seeks = Rc::new(Cell::new(0));
    let mut reference = TwoBitReference::new(TracedReader {
        inner: Cursor::new(raw),
        bytes_read: Rc::clone(&bytes_read),
        seeks: Rc::clone(&seeks),
    })
    .unwrap();
    bytes_read.set(0);
    seeks.set(0);
    let observed = reference.partition_base_counts::<4>("x", &runs).unwrap();
    assert_eq!(observed, expected);
    assert_eq!(bytes_read.get(), packed_bytes);
    assert_eq!(seeks.get(), 1);

    let mut combined = BaseCounts::default();
    for counts in observed {
        combined += counts;
    }
    assert_eq!(combined, reference.base_counts("x").unwrap());
}

#[test]
fn partition_validation_precedes_io_and_accepts_empty_contigs() {
    let bytes_read = Rc::new(Cell::new(0));
    let seeks = Rc::new(Cell::new(0));
    let mut reference = TwoBitReference::new(TracedReader {
        inner: Cursor::new(fixture(FixtureOrder::Little)),
        bytes_read: Rc::clone(&bytes_read),
        seeks: Rc::clone(&seeks),
    })
    .unwrap();
    bytes_read.set(0);
    seeks.set(0);
    let invalid: &[&[(u64, usize)]] = &[
        &[],
        &[(0, 0), (9, 0)],
        &[(5, 0), (5, 1), (9, 0)],
        &[(5, 0), (4, 1), (9, 0)],
        &[(8, 0)],
        &[(10, 0)],
        &[(1, 0), (9, 2)],
    ];
    for runs in invalid {
        assert!(reference.partition_base_counts::<2>("a", runs).is_err());
    }
    let absent = reference.partition_base_counts::<2>("absent", &[]);
    assert!(absent.is_err());
    let zero_width = reference.partition_base_counts::<2>("z", &[(0, 0)]);
    assert!(zero_width.is_err());
    let empty = reference.partition_base_counts::<2>("z", &[]).unwrap();
    assert_eq!(empty, [BaseCounts::default(); 2]);
    assert_eq!(bytes_read.get(), 0);
    assert_eq!(seeks.get(), 0);
}
