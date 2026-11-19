#!/usr/bin/env python3
"""Report public metadata for the OL27 sequence design and its IGVF reference."""
import hashlib
import json
import urllib.error
import urllib.request

ENCODE_ID = "ENCFF715XTT"
IGVF_ID = "IGVFFI9985FTZH"
ENCODE_URL = f"https://www.encodeproject.org/files/{ENCODE_ID}/?format=json"
IGVF_URL = f"https://api.data.igvf.org/tabular-files/{IGVF_ID}/?format=json"
CATALOG_URL = f"https://api.catalogkg.igvf.org/files-filesets?file_fileset_id={IGVF_ID}"
MAX_RESPONSE = 1024 * 1024
FIELDS = (
    "accession", "status", "content_type", "file_format", "file_format_type",
    "file_size", "md5sum", "content_md5sum", "sha256sum", "assembly",
    "output_type", "derived_from", "file_format_specifications",
    "data_use", "license", "controlled_access",
)


def identifiers(value):
    if isinstance(value, dict):
        return value.get("@id")
    if isinstance(value, list):
        return [identifiers(item) for item in value]
    return value


def public_urls(obj, fields):
    # Do not copy query parameters from URLs returned by a portal.
    return {key: obj[key].split("?", 1)[0] for key in fields
            if isinstance(obj.get(key), str)}


def project_object(obj, accession):
    if not isinstance(obj, dict) or obj.get("accession") != accession:
        raise ValueError("metadata accession disagrees")
    result = {key: identifiers(obj[key]) for key in FIELDS if key in obj}
    result.update(public_urls(obj, ("href", "s3_uri")))
    return result


def project_catalog(data):
    rows = data if isinstance(data, list) else [data]
    matches = [row for row in rows if isinstance(row, dict)
               and str(row.get("_id", "")).rsplit("/", 1)[-1] == IGVF_ID]
    if not matches:
        raise ValueError("catalog response lacks the expected file set")
    return [{"_id": row["_id"], "file_set_id": row["file_set_id"],
             **public_urls(row, ("source_url", "download_link"))}
            for row in matches]


def inspect(url, project):
    receipt = {"url": url}
    try:
        request = urllib.request.Request(
            url, headers={"Accept": "application/json", "User-Agent": "Helix-source-audit/1.0"}
        )
        with urllib.request.urlopen(request, timeout=20) as response:
            receipt["http_status"] = response.status
            if response.headers.get_content_type() not in ("application/json", "application/ld+json"):
                raise ValueError("expected JSON response")
            raw = response.read(MAX_RESPONSE + 1)
        if len(raw) > MAX_RESPONSE:
            raise ValueError("metadata response exceeds the read limit")
        receipt["response_sha256"] = hashlib.sha256(raw).hexdigest()
        receipt["metadata"] = project(json.loads(raw))
    except urllib.error.HTTPError as error:
        receipt["http_status"] = error.code
        receipt["failure"] = "HTTP response did not provide metadata"
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        receipt["failure"] = type(error).__name__
    except ValueError as error:
        receipt["failure"] = str(error)
    return receipt


def main():
    result = {
        "encode_design": inspect(ENCODE_URL, lambda obj: project_object(obj, ENCODE_ID)),
        "igvf_reference": inspect(IGVF_URL, lambda obj: project_object(obj, IGVF_ID)),
        "igvf_catalog": inspect(CATALOG_URL, project_catalog),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
