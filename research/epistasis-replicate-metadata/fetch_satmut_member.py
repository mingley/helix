#!/usr/bin/env python3
"""Fetch only the selected Zenodo design member by HTTP range, not outcomes."""
import argparse
import hashlib
import io
from pathlib import Path
import urllib.request
import zipfile

URL = 'https://zenodo.org/api/records/15297965/files/data_satmut.zip/content'
ARCHIVE_BYTES = 93320206
MEMBER = 'data/satmut/satmut_sequences.txt.gz'
MEMBER_BYTES = 3455761
MEMBER_COMPRESSED_BYTES = 2507667
MEMBER_CRC32 = 2316519572
MEMBER_SHA256 = 'dd0080c8d451710e766e6b2dadb35e952da680183d47eb5775adc7fbaf91dea0'


class RangeReader(io.RawIOBase):
    def __init__(self):
        self.position = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        if whence == 0:
            new = offset
        elif whence == 1:
            new = self.position + offset
        elif whence == 2:
            new = ARCHIVE_BYTES + offset
        else:
            raise ValueError('invalid seek mode')
        if not 0 <= new <= ARCHIVE_BYTES:
            raise ValueError('seek outside archive')
        self.position = new
        return new

    def read(self, count=-1):
        if count < 0:
            count = ARCHIVE_BYTES - self.position
        count = min(count, ARCHIVE_BYTES - self.position)
        if count <= 0:
            return b''
        end = self.position + count - 1
        request = urllib.request.Request(URL, headers={
            'Range': f'bytes={self.position}-{end}',
            'User-Agent': 'helix-research-metadata/1.0',
        })
        with urllib.request.urlopen(request, timeout=60) as response:
            expected = f'bytes {self.position}-{end}/{ARCHIVE_BYTES}'
            if response.status != 206 or response.headers.get('Content-Range') != expected:
                raise ValueError('server did not honor the requested ZIP range')
            result = response.read(count + 1)
        if len(result) != count:
            raise ValueError('ZIP range has unexpected length')
        self.position += count
        return result


def fetch():
    with zipfile.ZipFile(RangeReader()) as archive:
        member = archive.getinfo(MEMBER)
        if (member.file_size, member.compress_size, member.CRC) != (
                MEMBER_BYTES, MEMBER_COMPRESSED_BYTES, MEMBER_CRC32):
            raise ValueError('unexpected ZIP member metadata')
        # zipfile verifies CRC while reading the selected member. It does not
        # extract or decode any other archive member.
        data = archive.read(member)
    if len(data) != MEMBER_BYTES or hashlib.sha256(data).hexdigest() != MEMBER_SHA256:
        raise ValueError('satmut design member SHA-256 or size mismatch')
    return data


def write_new(path, data):
    with path.open('xb') as dest:
        dest.write(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='new design-member file')
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('output already exists')
    data = fetch()
    write_new(args.out, data)
    print(f'{MEMBER} bytes={len(data)} sha256={MEMBER_SHA256}')


if __name__ == '__main__':
    main()
