#!/usr/bin/env python3
"""Download species-wide paired-end and long-read RNA data from NCBI SRA.

Based on the standalone CLI conventions of 01_download_asmProt.py.
Requires Python >=3.9; downloads additionally require SRA Toolkit >=3.4.0.
See workflows/orthogroup_pangenome/README.md for usage and limitations.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import gzip
import hashlib
import http.client
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RUN_RE = re.compile(r"[SED]RR[0-9]+\Z")
READ_TYPES = frozenset(("paired_short", "pacbio", "nanopore"))
SRA_TABLE_FIELDS = (
    "Run", "Experiment", "SRAStudy", "BioProject", "Sample", "BioSample",
    "SampleName", "ScientificName", "TaxID", "read_type", "LibraryName",
    "LibraryStrategy", "LibrarySource", "LibrarySelection", "LibraryLayout",
    "Platform", "Model", "spots", "bases", "spots_with_mates", "avgLength",
    "size_MB", "ReleaseDate", "LoadDate", "CenterName", "Submission", "Consent",
    "selected", "selection_reason", "metadata_warning",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_json(path: Path, value: dict) -> None:
    """Publish small metadata atomically on the destination filesystem."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


class NcbiClient:
    def __init__(self, email=None, api_key=None, retries=3, timeout=120, interval=None):
        self.email = email
        self.api_key = api_key
        self.retries = retries
        self.timeout = timeout
        self.interval = interval if interval is not None else (0.11 if api_key else 0.35)
        self.last_request = 0.0

    def get(self, endpoint: str, **params) -> str:
        """POST E-utilities parameters; retry transient errors, respect rate limits.

        POST avoids URL-length limits for ID batches and keeps API keys out of
        URLs and exception messages. No credentials are written to provenance.
        """
        data = {**params, "tool": "species_rnaseq_downloader"}
        if self.email:
            data["email"] = self.email
        if self.api_key:
            data["api_key"] = self.api_key
        request = urllib.request.Request(
            EUTILS + endpoint, data=urllib.parse.urlencode(data).encode("utf-8"),
            headers={"User-Agent": "species-rnaseq-downloader/1.0"},
        )
        for attempt in range(self.retries):
            time.sleep(max(0, self.interval - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read().decode("utf-8-sig")
            except urllib.error.HTTPError as exc:
                retry = exc.code == 429 or 500 <= exc.code < 600
                message = f"NCBI {endpoint}: HTTP {exc.code}"
            except (urllib.error.URLError, TimeoutError, OSError, http.client.HTTPException) as exc:
                retry = True
                message = f"NCBI {endpoint}: {type(exc).__name__}"
            if not retry or attempt + 1 == self.retries:
                raise RuntimeError(message) from None
            time.sleep(min(2 ** attempt, 30))
        raise RuntimeError("No HTTP attempts configured")


def search_result(text: str) -> dict:
    try:
        result = json.loads(text)
        if "error" in result:
            raise ValueError(str(result["error"]))
        result = result["esearchresult"]
        if result.get("errorlist"):
            raise ValueError(str(result["errorlist"]))
        int(result["count"])
        return result
    except (ValueError, KeyError, TypeError) as exc:
        raise RuntimeError(f"Invalid NCBI search response: {exc}") from exc


def resolve_taxon(client, species=None, taxid=None) -> dict:
    if bool(species) == bool(taxid):
        raise RuntimeError("Supply a species name or --taxid, exclusively")
    if species:
        if any(ch in species for ch in '"[]\n\r'):
            raise RuntimeError("Species must be a name, not an Entrez query")
        result = search_result(client.get("esearch.fcgi", db="taxonomy",
                               term=f'"{species}"[Scientific Name]', retmode="json", retmax=2))
        count = int(result["count"])
        if count == 0:
            raise RuntimeError(f"No NCBI taxon found for {species!r}; try --taxid")
        if count != 1 or len(result.get("idlist", [])) != 1:
            raise RuntimeError(f"Ambiguous species name {species!r}; use an explicit --taxid")
        taxid = result["idlist"][0]
    taxid = str(taxid)
    if not taxid.isdigit() or int(taxid) <= 0:
        raise RuntimeError("TaxID must be a positive integer")
    try:
        xml = ET.fromstring(client.get("efetch.fcgi", db="taxonomy", id=taxid, retmode="xml"))
    except ET.ParseError as exc:
        raise RuntimeError("NCBI returned invalid taxonomy XML") from exc
    taxa = xml.findall("./Taxon")
    if len(taxa) != 1 or not taxa[0].findtext("ScientificName"):
        raise RuntimeError(f"No unambiguous taxonomy record for TaxID {taxid}")
    node = taxa[0]
    return {"taxid": node.findtext("TaxId"), "name": node.findtext("ScientificName"),
            "rank": node.findtext("Rank", "unknown")}


def build_query(taxid: str, exact_taxon=False) -> str:
    taxid = str(taxid)
    if not taxid.isdigit():
        raise ValueError("Invalid TaxID")
    organism = f"txid{taxid}[Organism:noexp]" if exact_taxon else f"txid{taxid}[Organism]"
    # Broader than eligibility: missing RNA source must not hide labelled RNA-Seq.
    # The RNA property also captures long reads classified outside RNA-Seq.
    return (
        f'{organism} AND public[Access] AND ("rna seq"[Strategy] OR '
        '((pacbio_smrt[Platform] OR oxford_nanopore[Platform]) AND '
        '(TRANSCRIPTOMIC[Source] OR biomol_rna[Properties])))'
    )


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def classify_run(row: dict) -> tuple[str, str]:
    consent = normalized(row.get("Consent", ""))
    if consent and consent != "public":
        return "", "excluded: non-public/controlled-access consent"
    platform = normalized(row.get("Platform", ""))
    strategy = normalized(row.get("LibraryStrategy", ""))
    source = normalized(row.get("LibrarySource", ""))
    layout = row.get("LibraryLayout", "").strip().upper()
    long_type = {"pacbiosmrt": "pacbio", "pacbio": "pacbio",
                 "oxfordnanopore": "nanopore"}.get(platform)
    if long_type and (strategy == "rnaseq" or source == "transcriptomic"):
        return long_type, "long-read RNA-Seq or TRANSCRIPTOMIC source"
    if strategy == "rnaseq" and layout == "PAIRED":
        return "paired_short", "paired-end RNA-Seq"
    if not strategy or not platform or not layout:
        return "", "excluded: missing strategy, platform, or layout metadata"
    return "", "excluded: outside paired RNA-Seq / long transcriptomic scope"


def write_manifest(root: Path, rows: list[dict]) -> None:
    core = ["Run", "Experiment", "read_type", "selected", "selection_reason", "metadata_warning"]
    fields = core + sorted({key for row in rows for key in row} - set(core))
    with (root / "runs.tsv.tmp").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (root / "runs.tsv.tmp").replace(root / "runs.tsv")
    (root / "selected_runs.txt").write_text(
        "".join(row["Run"] + "\n" for row in rows if row["selected"] == "yes"), encoding="utf-8")


def write_sra_run_table(root: Path, rows: list[dict]) -> Path:
    """Export selected RunInfo records, keeping every supplied metadata field.

    This is a RunInfo-derived table, not NCBI's browser SraRunTable.txt export;
    custom sample attributes absent from RunInfo cannot be reconstructed here.
    """
    selected = [row for row in rows if row.get("selected") == "yes"]
    fields = list(SRA_TABLE_FIELDS) + sorted({key for row in selected for key in row} - set(SRA_TABLE_FIELDS))
    root.mkdir(parents=True, exist_ok=True)
    path = root / "sra_run_table.tsv"
    temporary = path.with_suffix(".tsv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    temporary.replace(path)
    return path


def print_sra_run_table(rows: list[dict], path: Path, limit=20) -> None:
    """Print a bounded preview; the TSV always contains every selected row."""
    if limit < 0:
        raise ValueError("List limit must be zero or positive")
    selected = [row for row in rows if row.get("selected") == "yes"]
    preview = selected if limit == 0 else selected[:limit]
    print(f"\nSRA selector-style table (RunInfo; showing {len(preview)} of {len(selected)} selected runs)")
    columns = ("Run", "BioProject", "BioSample", "SampleName", "ScientificName", "read_type",
               "LibraryStrategy", "LibraryLayout", "Platform", "Model", "spots", "bases", "size_MB")
    if preview:
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in preview:
            # Keep raw TSV metadata untouched, but don't let embedded tabs or
            # newlines in sample descriptions split the terminal preview rows.
            values = []
            for field in columns:
                value = row.get(field)
                value = "" if value is None else str(value)
                values.append(re.sub(r"[\x00-\x1f\x7f]+", " ", value).strip() or "-")
            writer.writerow(values)
    else:
        print("No selected runs.")
    print("RunInfo-derived; not an exact copy of NCBI's browser SraRunTable.txt. Only attributes in RunInfo are shown.")
    print("spots are sequencing spots, not individual reads; size_MB is reported archive size, not FASTQ size.")
    print(f"Full selected-run metadata ({len(selected)} rows, all available columns): {path}")
    if len(preview) < len(selected):
        print("Use --list-limit 0 to print every selected run; the TSV is never preview-limited.")


def discover_runs(client, taxon: dict, root: Path, exact_taxon=False, page_size=500, read_types=None) -> tuple[list[dict], dict]:
    if not 1 <= page_size <= 10000:
        raise ValueError("Page size must be 1..10000")
    enabled = READ_TYPES if read_types is None else frozenset(read_types)
    if not enabled or not enabled.issubset(READ_TYPES):
        raise ValueError("Select at least one valid read type: paired_short, pacbio, nanopore")
    root.mkdir(parents=True, exist_ok=True)
    raw_dir = root / "metadata" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    raw_dir.mkdir(parents=True)
    query = build_query(taxon["taxid"], exact_taxon)
    info = {"taxon": taxon, "query": query, "retrieved_utc": utc_now(), "complete": False,
            "raw_metadata_dir": str(raw_dir), "experiment_count": 0, "run_count": 0,
            "read_types": sorted(enabled)}
    write_json(root / "discovery.json", info)
    try:
        raw_search = client.get("esearch.fcgi", db="sra", term=query, retmax=0,
                                usehistory="y", retmode="json")
        (raw_dir / "search.json").write_text(raw_search, encoding="utf-8")
        search = search_result(raw_search)
        count = int(search["count"])
        info["experiment_count"] = count
        if count and (not search.get("webenv") or not search.get("querykey")):
            raise RuntimeError("NCBI search returned no history snapshot")
        info["history"] = {"WebEnv": search.get("webenv"), "query_key": search.get("querykey")}
        seen_ids, seen_experiments, runs = set(), set(), {}
        for start in range(0, count, page_size):
            print(f"Metadata: experiments {start + 1}–{min(start + page_size, count)} / {count}", flush=True)
            raw_page = client.get("esearch.fcgi", db="sra", term=f'#{search["querykey"]}',
                                  WebEnv=search["webenv"], usehistory="y", retstart=start, retmax=page_size, retmode="json")
            (raw_dir / f"ids_{start:09d}.json").write_text(raw_page, encoding="utf-8")
            page = search_result(raw_page)
            ids = page.get("idlist", [])
            if (int(page["count"]) != count or len(ids) != min(page_size, count - start)
                    or len(set(ids)) != len(ids) or seen_ids.intersection(ids)):
                raise RuntimeError("Incomplete or duplicate history page; search coverage cannot be verified")
            seen_ids.update(ids)
            raw = client.get("efetch.fcgi", db="sra", id=",".join(ids), rettype="runinfo", retmode="text")
            (raw_dir / f"runinfo_{start:09d}.csv").write_text(raw, encoding="utf-8")
            reader = csv.DictReader(io.StringIO(raw.lstrip()))
            if not reader.fieldnames or not {"Run", "Experiment"}.issubset(reader.fieldnames):
                raise RuntimeError("NCBI returned invalid RunInfo CSV")
            experiments = set()
            for item in reader:
                if None in item or any(value is None for value in item.values()):
                    raise RuntimeError("Truncated or malformed RunInfo row")
                item = {k: v.strip() for k, v in item.items()}
                run, experiment = item["Run"], item["Experiment"]
                if not RUN_RE.fullmatch(run) or not experiment:
                    raise RuntimeError("RunInfo contains an invalid run/experiment accession")
                experiments.add(experiment)
                if run in runs and runs[run] != item:
                    raise RuntimeError(f"Conflicting metadata for run {run}")
                runs[run] = item
            if len(experiments) != len(ids) or experiments.intersection(seen_experiments):
                raise RuntimeError("Experiment coverage mismatch in RunInfo; discovery is incomplete")
            seen_experiments.update(experiments)
        rows = []
        for run in sorted(runs):
            item = runs[run]
            bucket, reason = classify_run(item)
            if bucket and bucket not in enabled:
                reason = f"excluded by requested read-type filter: {reason}"
            warnings = []
            if normalized(item.get("LibraryStrategy", "")) == "rnaseq" and normalized(item.get("LibrarySource", "")) != "transcriptomic":
                warnings.append("RNA-Seq label with missing/non-TRANSCRIPTOMIC source")
            if any(not item.get(k) for k in ("Platform", "LibraryLayout", "LibraryStrategy", "LibrarySource")):
                warnings.append("missing classification metadata")
            rows.append({**item, "read_type": bucket, "selected": "yes" if bucket in enabled else "no",
                         "selection_reason": reason, "metadata_warning": "; ".join(warnings)})
        info.update(complete=True, run_count=len(rows), selected_run_count=sum(row["selected"] == "yes" for row in rows),
                    retrieved_experiment_count=len(seen_experiments))
        write_manifest(root, rows)
        write_json(root / "discovery.json", info)
        return rows, info
    except Exception as exc:
        info["error"] = str(exc)
        write_json(root / "discovery.json", info)
        raise


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_completion(run_dir: Path, settings: dict) -> bool:
    """A completion marker is trusted only with identical settings and bytes."""
    try:
        record = json.loads((run_dir / "complete.json").read_text(encoding="utf-8"))
        if (not isinstance(record, dict) or record.get("Run") != run_dir.name
                or record.get("status") != "complete" or record.get("settings") != settings
                or record.get("archive_format") != "SRA Normalized"
                or not isinstance(record.get("files"), list) or not record["files"]):
            return False
        names = set()
        for item in record["files"]:
            if not isinstance(item, dict):
                return False
            name = item["name"]
            if (not isinstance(name, str) or Path(name).name != name or name in names
                    or not re.fullmatch(re.escape(run_dir.name) + r"(?:_[0-9]+)?\.fastq\.gz", name)):
                return False
            names.add(name)
            path = run_dir / name
            if not path.is_file() or path.stat().st_size != item["size"] or sha256_file(path) != item["sha256"]:
                return False
        return names == {path.name for path in run_dir.glob("*.fastq.gz")}
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_command(command: list[str], log_path: Path, env: dict, retries=1) -> str:
    """Keep potentially long-running Toolkit output on disk, not in memory."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        with log_path.open("ab") as log:
            log.write((f"\n[{utc_now()}] attempt {attempt + 1}: "
                       + repr([str(arg) for arg in command]) + "\n").encode("utf-8"))
            log.flush()
            start = log.tell()
            try:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=False)
            except OSError as exc:
                raise RuntimeError(f"Cannot execute {command[0]}: {exc}") from exc
        with log_path.open("rb") as log:
            log.seek(max(start, log_path.stat().st_size - 1024 * 1024))
            output = log.read().decode("utf-8", errors="replace")
        if result.returncode == 0:
            return output
        if attempt + 1 < retries:
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"{command[0]} exited {result.returncode}; see {log_path}")


def download_run(row: dict, root: Path, cache_dir: Path, tmp_dir: Path,
                 threads=4, retries=3, pacbio_table="auto", include_technical=False, tool_versions=None) -> dict:
    """Resume an archive, validate it, then publish one checksummed FASTQ set."""
    run = row["Run"]
    if not RUN_RE.fullmatch(run):
        raise RuntimeError(f"Invalid run accession: {run!r}")
    if pacbio_table not in ("auto", "SEQUENCE") or threads < 1 or retries < 1:
        raise ValueError("Invalid conversion options")
    settings = {"pacbio_table": pacbio_table, "include_technical": include_technical}
    fastq_root = root / "fastq"
    fastq_root.mkdir(parents=True, exist_ok=True)
    run_dir = fastq_root / run
    if verify_completion(run_dir, settings):
        record = json.loads((run_dir / "complete.json").read_text(encoding="utf-8"))
        return {**record, "resumed": True}
    cache = cache_dir / run
    cache.mkdir(parents=True, exist_ok=True)
    archive = cache / f"{run}.sra"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    log_path = root / "logs" / f"{run}.log"
    env = dict(os.environ)
    # Process-local preference; never edit ~/.ncbi/user-settings.mkfg.
    env["NCBI_VDB_QUALITY"] = "R"
    print(f"[{run}] prefetch (normalized SRA; --max-size u), then validate", flush=True)
    run_command(["prefetch", run, "--type", "sra", "--max-size", "u", "-O", str(cache)],
                log_path, env, retries=retries)
    if not archive.is_file() or archive.stat().st_size == 0:
        raise RuntimeError(f"No normalized archive at {archive}; SRA Lite is not accepted")
    run_command(["vdb-validate", str(archive)], log_path, env)
    with tempfile.TemporaryDirectory(prefix=f".{run}-", dir=fastq_root) as stage_name, \
            tempfile.TemporaryDirectory(prefix=f"{run}-", dir=tmp_dir) as scratch_name:
        stage = Path(stage_name)
        command = ["fasterq-dump", str(cache), "--outdir", str(stage), "--temp", scratch_name,
                   "--threads", str(threads), "--details"]
        command += ["--split-files", "--include-technical"] if include_technical else ["--split-3"]
        is_pacbio = normalized(row.get("Platform", "")) in ("pacbio", "pacbiosmrt")
        if is_pacbio and pacbio_table == "SEQUENCE":
            command += ["--table", "SEQUENCE"]
        print(f"[{run}] convert FASTQ ({threads} threads), then gzip", flush=True)
        details = run_command(command, log_path, env)
        files = []
        for fastq in sorted(stage.glob("*.fastq")):
            if not re.fullmatch(re.escape(run) + r"(?:_[0-9]+)?\.fastq", fastq.name):
                raise RuntimeError(f"Unexpected FASTQ output name: {fastq.name}")
            if not fastq.stat().st_size:
                continue
            compressed = fastq.with_suffix(".fastq.gz")
            with fastq.open("rb") as source, gzip.open(compressed, "wb", compresslevel=1) as dest:
                shutil.copyfileobj(source, dest, 1024 * 1024)
            files.append({"name": compressed.name, "size": compressed.stat().st_size,
                          "sha256": sha256_file(compressed)})
            fastq.unlink()  # Only a temporary conversion output owned by this invocation.
        if not files:
            raise RuntimeError(f"{run}: conversion produced no nonempty FASTQ files")
        record = {"Run": run, "status": "complete", "completed_utc": utc_now(),
                  "settings": settings, "archive_format": "SRA Normalized",
                  "archive_path": str(archive), "archive_size": archive.stat().st_size,
                  "quality_preference": "R", "source_metadata": row, "files": files,
                  "tool_versions": tool_versions or {},
                  "conversion_details": details,
                  "table_policy": "PacBio SEQUENCE override" if is_pacbio and pacbio_table == "SEQUENCE"
                  else ("PacBio default: consensus when available, otherwise sequence" if is_pacbio
                        else "Toolkit default (see conversion details)")}
        write_json(stage / "complete.json", record)
        # Retain old output on a changed/corrupt rerun. Do not mingle stale files
        # with a new output set, or delete files supplied by the user.
        previous = None
        if run_dir.exists():
            previous = fastq_root / f"{run}.previous-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}"
            run_dir.rename(previous)
        try:
            stage.rename(run_dir)
        except OSError:
            if previous is not None:
                previous.rename(run_dir)
            raise
        return record


def check_toolkit() -> dict:
    versions = {}
    for tool in ("prefetch", "vdb-validate", "fasterq-dump"):
        if not shutil.which(tool):
            raise RuntimeError(f"SRA Toolkit executable {tool!r} not found. Install sra-tools >=3.4.0: https://github.com/ncbi/sra-tools")
        result = subprocess.run([tool, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, check=False)
        match = re.search(r"\b([0-9]+)\.([0-9]+)\.([0-9]+)\b", result.stdout)
        if result.returncode or not match or tuple(map(int, match.groups())) < (3, 4, 0):
            raise RuntimeError(f"{tool}: SRA Toolkit >=3.4.0 is required for long reads; found {result.stdout.strip()!r}")
        versions[tool] = match.group(0)
    return versions


@contextmanager
def project_lock(root: Path):
    """A kernel-released lock avoids stale locks after interrupted cluster jobs."""
    # Intended for Linux/HPC and macOS, like the assembly downloader's wrapper.
    import fcntl
    with (root / ".rnaseq.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another downloader is using {root}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def positive_int(value: str) -> int:
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download all discoverable paired-end RNA-Seq and PacBio/ONT RNA runs for a species from NCBI.",
        epilog="Includes strains/subspecies by default. --list-only transfers metadata only. See workflows/orthogroup_pangenome/README.md.",
    )
    parser.add_argument("species", nargs="?", help="Scientific species name, quoted if it contains spaces")
    parser.add_argument("--taxid", type=positive_int, help="NCBI TaxID instead of a species name")
    parser.add_argument("--out-dir", type=Path, default=Path("data/rnaseq"), metavar="DIR", help="Output parent (default: data/rnaseq)")
    parser.add_argument("--dl-dir", type=Path, metavar="DIR", help="Archive cache parent (default: species-output/_sra)")
    parser.add_argument("--tmp-dir", type=Path, metavar="DIR", help="Conversion scratch parent (default: species-output/_tmp)")
    parser.add_argument("--threads", type=positive_int, default=4, help="Threads per conversion; runs processed sequentially (default: 4)")
    parser.add_argument("--retries", type=positive_int, default=3, help="Maximum HTTP/prefetch attempts (default: 3)")
    parser.add_argument("--timeout", type=positive_int, default=120, help="HTTP timeout seconds (default: 120)")
    parser.add_argument("--list-only", action="store_true", help="Preview RunInfo metadata and export sra_run_table.tsv without downloading reads or needing SRA Toolkit. RunInfo-derived; not an exact copy of NCBI's browser SraRunTable.txt")
    parser.add_argument("--list-limit", type=int, default=20, metavar="N", help="Maximum selected rows printed by --list-only (default: 20; 0 prints all); TSV always includes all selected runs")
    parser.add_argument("--paired-end", action="store_true", help="Select paired-end short RNA-Seq; combine with --long-reads for both (default: both when neither flag is given)")
    parser.add_argument("--long-reads", action="store_true", help="Select PacBio/ONT RNA reads regardless of layout; combine with --paired-end for both")
    parser.add_argument("--exact-taxon", action="store_true", help="Only the exact TaxID, excluding descendant strain/subspecies nodes")
    parser.add_argument("--pacbio-table", choices=("auto", "SEQUENCE"), default="auto", help="PacBio table: consensus when available (auto), or SEQUENCE; does not override ONT")
    parser.add_argument("--include-technical", action="store_true", help="Include barcode/technical reads using split-files; disables split-3 orphan separation")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"), help="Contact email for NCBI (or NCBI_EMAIL)")
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.species) == bool(args.taxid):
        parser.error("Supply a species name or --taxid, exclusively")
    if args.list_limit < 0:
        parser.error("--list-limit must be zero or positive")
    read_types = set()
    if args.paired_end:
        read_types.add("paired_short")
    if args.long_reads:
        read_types.update(("pacbio", "nanopore"))
    if not read_types:
        read_types.update(READ_TYPES)
    client = NcbiClient(email=args.email, api_key=os.environ.get("NCBI_API_KEY"),
                        retries=args.retries, timeout=args.timeout)
    try:
        taxon = resolve_taxon(client, species=args.species, taxid=args.taxid)
        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", taxon["name"]).strip("._-") or "taxon"
        root = args.out_dir.expanduser().resolve() / f'{safe_name}_taxid{taxon["taxid"]}'
        root.mkdir(parents=True, exist_ok=True)
        with project_lock(root):
            print(f'Taxon: {taxon["name"]} (TaxID {taxon["taxid"]}, rank {taxon["rank"]})\nOutput: {root}', flush=True)
            rows, info = discover_runs(client, taxon, root, exact_taxon=args.exact_taxon, read_types=read_types)
            if not info.get("complete"):
                raise RuntimeError("Discovery did not complete; downloads will not start")
            selected = [row for row in rows if row["selected"] == "yes"]
            counts = Counter(row["read_type"] for row in selected)
            print(f"Discovered {len(rows)} runs; selected {len(selected)}: "
                  + ", ".join(f"{kind}={number}" for kind, number in sorted(counts.items())))
            known_sizes = []
            for row in selected:
                try:
                    known_sizes.append(float(row["size_MB"]))
                except (ValueError, KeyError):
                    pass
            if known_sizes:
                print(f"RunInfo archive-size estimate: {sum(known_sizes):,.1f} MB for {len(known_sizes)}/{len(selected)} runs; normalized SRA may differ, and FASTQ/scratch need additional space.")
            if args.list_only:
                table_path = write_sra_run_table(root, rows)
                print_sra_run_table(rows, table_path, limit=args.list_limit)
            if args.list_only or not selected:
                print("Metadata saved; no sequencing reads downloaded.")
                return 0
            versions = check_toolkit()
            cache = args.dl_dir.expanduser().resolve() if args.dl_dir else root / "_sra"
            scratch = args.tmp_dir.expanduser().resolve() if args.tmp_dir else root / "_tmp"
            execution = {"started_utc": utc_now(), "taxon": taxon, "tool_versions": versions,
                         "read_types": sorted(read_types),
                         "python": sys.version.split()[0], "cache_dir": str(cache), "tmp_dir": str(scratch),
                         "settings": {"pacbio_table": args.pacbio_table, "include_technical": args.include_technical},
                         "threads": args.threads, "quality_preference": "R", "prefetch_max_size": "u",
                         "runs": {row["Run"]: {"status": "pending"} for row in selected}}
            status_path = root / "download_status.json"
            write_json(status_path, execution)
            for index, row in enumerate(selected, 1):
                run = row["Run"]
                print(f"[{index}/{len(selected)}] {run} ({row['read_type']})", flush=True)
                try:
                    result = download_run(row, root, cache, scratch, threads=args.threads, retries=args.retries,
                                          pacbio_table=args.pacbio_table, include_technical=args.include_technical,
                                          tool_versions=versions)
                    execution["runs"][run] = {"status": result["status"], "resumed": result.get("resumed", False),
                                               "files": result["files"], "tool_versions": result.get("tool_versions", {})}
                except (RuntimeError, OSError) as exc:
                    execution["runs"][run] = {"status": "failed", "error": str(exc)}
                    print(f"[{run}] FAILED: {exc}", file=sys.stderr, flush=True)
                write_json(status_path, execution)
            failures = sum(record["status"] != "complete" for record in execution["runs"].values())
            execution.update(finished_utc=utc_now(), complete=failures == 0)
            write_json(status_path, execution)
            print(f"Finished: {len(selected) - failures} complete, {failures} failed. Status: {status_path}")
            return 1 if failures else 0
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted. Rerun the same command to resume validated outputs and archive downloads.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
