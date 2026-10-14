"""
Download Ensembl reference genomes and register them in genomes.yaml.

Supported sources (HTTPS views of the Ensembl FTP sites):

* Ensembl (vertebrates plus the model organisms it hosts, e.g. yeast, worm, fly):
  ``https://ftp.ensembl.org/pub/release-<N>/``
* Ensembl Genomes plants, fungi, metazoa and protists, including species kept in
  ``*_collection`` directories:
  ``https://ftp.ensemblgenomes.ebi.ac.uk/pub/<division>/release-<N>/``

Ensembl Bacteria is not supported.

Downloads are verified against the ``CHECKSUMS`` file of their directory.
Ensembl writes it with the BSD ``sum`` algorithm (16-bit rotating checksum
and the size in 1 KiB blocks), e.g. ``22450 861294 <file name>``.
"""

import datetime
import difflib
import gzip
import http.client
import json
import os
import re
import shutil
import subprocess
import time
import typing as T
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import yaml

from uprobe import __version__
from .utils import get_logger

log = get_logger(__name__)

ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
ENSEMBL_GENOMES_FTP = "https://ftp.ensemblgenomes.ebi.ac.uk/pub"
ENSEMBL_REST = "https://rest.ensembl.org"
# CLI division name -> Ensembl division id (as used in species_<id>.txt).
DIVISIONS = {
    "vertebrates": "EnsemblVertebrates",
    "plants": "EnsemblPlants",
    "fungi": "EnsemblFungi",
    "metazoa": "EnsemblMetazoa",
    "protists": "EnsemblProtists",
}
MANIFEST_NAME = "genome_source.json"
CHECKSUM_ALGORITHM = "BSD sum (16-bit checksum, 1 KiB blocks)"
# Files U-Probe derives from the FASTA/GTF next to them (see uprobe.core.tools).
_DERIVED = ("bowtie2_genome", "blast_genome", "mmseqs_genome",
            "bowtie2_transcript", "blast_transcript", "jf_genome", "transcript.fa")
# Rough "gz + decompressed" size factors, only used for a disk-space warning.
_SPACE_FACTOR = {"fasta": 5, "gtf": 12}
_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
_CHUNK = 1 << 20
_TIMEOUT = 60
_ATTEMPTS = 3
# Below this size the pure-Python checksum is fast enough (and the system `sum`
# variant cannot be told apart by its block count).
_PY_SUM_LIMIT = 1 << 20


class GenomeFetchError(RuntimeError):
    """A genome could not be resolved, downloaded or verified."""


@dataclass
class RemoteFile:
    kind: str                      # "fasta" or "gtf"
    url: str
    name: str                      # remote (.gz) file name
    checksum: T.Tuple[int, int]    # BSD sum: (checksum, 1 KiB blocks)

    @property
    def local_name(self) -> str:
        return self.name[:-3] if self.name.endswith(".gz") else self.name


@dataclass
class FetchPlan:
    division: str
    release: int
    species: str
    display_name: str
    assembly: str
    sequence: str                  # "primary_assembly" or "toplevel"
    fasta: RemoteFile
    gtf: RemoteFile

    def to_dict(self) -> dict:
        return {
            "division": self.division,
            "release": self.release,
            "species": self.species,
            "display_name": self.display_name,
            "assembly": self.assembly,
            "sequence": self.sequence,
            "files": {
                f.kind: {"url": f.url, "checksum": _fmt_checksum(f.checksum),
                         "compressed_mib": round(f.checksum[1] / 1024, 1)}
                for f in (self.fasta, self.gtf)
            },
        }


# ---------------------------------------------------------------- HTTP helpers

def _open(url: str, offset: int = 0):
    headers = {"User-Agent": f"uprobe/{__version__}"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    return urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=_TIMEOUT)


def _get_text(url: str, missing_ok: bool = False, attempts: int = _ATTEMPTS) -> T.Optional[str]:
    """GET a small text file; return None for HTTP 404 when ``missing_ok``."""
    for attempt in range(1, attempts + 1):
        try:
            with _open(url) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 404 and missing_ok:
                return None
            raise GenomeFetchError(f"HTTP {e.code} for {url}") from e
        except (OSError, http.client.HTTPException) as e:
            if attempt == attempts:
                raise GenomeFetchError(
                    f"Cannot reach {url}: {e}. Check the network connection, or use --mirror."
                ) from e
            time.sleep(attempt)


def _base_url(division: str, mirror: T.Optional[str] = None) -> str:
    if division not in DIVISIONS:
        raise GenomeFetchError(
            f"Unsupported division '{division}'. Choose from: {', '.join(DIVISIONS)} "
            "(Ensembl Bacteria is not supported)."
        )
    return (mirror or (ENSEMBL_FTP if division == "vertebrates" else ENSEMBL_GENOMES_FTP)).rstrip("/")


def release_root(division: str, release: int, mirror: T.Optional[str] = None) -> str:
    base = _base_url(division, mirror)
    if division == "vertebrates":
        return f"{base}/release-{release}"
    return f"{base}/{division}/release-{release}"


# ------------------------------------------------------ release and species

def resolve_release(release: T.Union[int, str, None] = "latest", division: str = "vertebrates",
                    mirror: T.Optional[str] = None, rest_url: T.Optional[str] = ENSEMBL_REST) -> int:
    """Return an explicit release number; ``latest`` is read from the FTP
    ``VERSION`` file, falling back to Ensembl REST."""
    if release is not None and str(release).strip().lower() not in ("latest", "current"):
        text = str(release).strip().lower().removeprefix("release-")
        if not text.isdigit() or int(text) <= 0:
            raise GenomeFetchError(f"Invalid release '{release}': use a release number (e.g. 116) or 'latest'.")
        return int(text)
    base = _base_url(division, mirror)
    try:
        text = _get_text(f"{base}/VERSION", missing_ok=True)
        if text and text.strip().isdigit():
            return int(text.strip())
    except GenomeFetchError as e:
        log.warning(f"Could not read {base}/VERSION: {e}")
    if rest_url:
        # Ensembl Genomes divisions share one release number (eg_version).
        endpoint = "info/data" if division == "vertebrates" else "info/eg_version"
        try:
            data = json.loads(_get_text(f"{rest_url.rstrip('/')}/{endpoint}?content-type=application/json",
                                        attempts=1))
            return int(max(data["releases"]) if "releases" in data else data["version"])
        except (GenomeFetchError, ValueError, KeyError, TypeError) as e:
            log.warning(f"Could not resolve the latest release from Ensembl REST: {e}")
    raise GenomeFetchError("Could not resolve the latest Ensembl release; pass an explicit --release number.")


def load_species_table(division: str, release: int, mirror: T.Optional[str] = None) -> T.List[dict]:
    """Rows of the release's ``species_<Division>.txt`` ([] if it is absent)."""
    url = f"{release_root(division, release, mirror)}/species_{DIVISIONS[division]}.txt"
    text = _get_text(url, missing_ok=True)
    if text is None:
        log.warning(f"Species list not found ({url}); using the species name as given.")
        return []
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    header = lines[0].lstrip("#").split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if not line.startswith("#")]


def _norm(text) -> str:
    return re.sub(r"[\s\-]+", "_", str(text).strip().lower())


def search_species(query: str, table: T.List[dict], limit: int = 10) -> T.List[dict]:
    """Species whose name, display name or taxonomy id matches ``query``,
    followed by close spellings."""
    q = _norm(query)
    hits = [r for r in table
            if q in r.get("species", "") or q in _norm(r.get("name", "")) or q == r.get("taxonomy_id")]
    hits.sort(key=lambda r: (q not in (r.get("species"), _norm(r.get("name", ""))),
                             not r.get("species", "").startswith(q), len(r.get("species", ""))))
    labels = {}
    for r in table:
        labels.setdefault(r.get("species", ""), r)
        labels.setdefault(_norm(r.get("name", "")), r)
    for label in difflib.get_close_matches(q, list(labels), n=limit, cutoff=0.6):
        if labels[label] not in hits:
            hits.append(labels[label])
    return hits[:limit]


def _rest_division_hint(q: str, division: str, rest_url: T.Optional[str]) -> T.Optional[str]:
    """Best-effort hint from Ensembl REST (it also knows aliases and other divisions)."""
    if not rest_url:
        return None
    try:
        data = json.loads(_get_text(
            f"{rest_url.rstrip('/')}/info/genomes/{urllib.parse.quote(q)}?content-type=application/json",
            attempts=1))
    except (GenomeFetchError, ValueError):
        return None
    if not isinstance(data, dict):
        return None
    found = {v: k for k, v in DIVISIONS.items()}.get(data.get("division"))
    name = data.get("name")
    if data.get("division") == "EnsemblBacteria":
        return "It is an Ensembl Bacteria genome, which is not supported."
    if found and found != division:
        return f"Ensembl REST lists '{name}' in {data['division']}: retry with --division {found} --species {name}."
    if found and name and name != q:
        return f"Ensembl REST resolves it to '{name}': retry with --species {name}."
    return None


def resolve_species(query: str, table: T.List[dict], division: str = "vertebrates",
                    release: T.Optional[int] = None, rest_url: T.Optional[str] = ENSEMBL_REST) -> dict:
    """Match a production name (homo_sapiens), display name (Human) or
    taxonomy id (9606) against the species table."""
    q = _norm(query)
    if not table:
        return {"species": q, "name": query}
    matches = [r for r in table if r.get("species") == q]
    if not matches:
        matches = [r for r in table if _norm(r.get("name", "")) == q]
    if not matches and q.isdigit():
        matches = [r for r in table if r.get("taxonomy_id") == q]
    if len(matches) == 1:
        return matches[0]
    where = f"Ensembl {division} release {release}"
    if matches:
        names = ", ".join(r["species"] for r in matches[:10])
        raise GenomeFetchError(f"'{query}' matches several species in {where}: {names}. Pass one of them to --species.")
    msg = [f"Species '{query}' not found in {where}."]
    suggestions = search_species(query, table, limit=5)
    if suggestions:
        msg.append("Did you mean: " + ", ".join(f"{r['species']} ({r.get('name')})" for r in suggestions) + "?")
    hint = _rest_division_hint(q, division, rest_url)
    if hint:
        msg.append(hint)
    msg.append(f"Search with `uprobe genome search <term> --division {division}`; "
               "non-vertebrates need --division plants|fungi|metazoa|protists.")
    raise GenomeFetchError(" ".join(msg))


def _species_path(row: dict) -> str:
    # Ensembl Genomes keeps many species in <division>_<clade>_collection/ directories.
    core_db = row.get("core_db", "")
    if "_collection_core_" in core_db:
        return f"{core_db.split('_core_')[0]}/{row['species']}"
    return row["species"]


# ------------------------------------------------------------- checksums

def parse_checksums(text: str) -> T.Dict[str, T.Tuple[int, int]]:
    """Parse a CHECKSUMS file into ``{file name: (checksum, blocks)}``."""
    result = {}
    for line in text.splitlines():
        parts = line.split(None, 2)
        if len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit():
            result[parts[2].strip()] = (int(parts[0]), int(parts[1]))
    return result


def _fmt_checksum(value: T.Tuple[int, int]) -> str:
    return f"{value[0]:05d} {value[1]}"


def _bsd_sum_python(path: Path) -> T.Tuple[int, int]:
    checksum = size = 0
    with open(path, "rb") as fh:
        while chunk := fh.read(_CHUNK):
            size += len(chunk)
            for byte in chunk:
                checksum = (((checksum >> 1) | ((checksum & 1) << 15)) + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def bsd_sum(path: T.Union[str, Path]) -> T.Tuple[int, int]:
    """BSD ``sum`` of a file: (16-bit checksum, size in 1 KiB blocks)."""
    path = Path(path)
    size = path.stat().st_size
    blocks = (size + 1023) // 1024
    exe = shutil.which("sum")
    if exe and size >= _PY_SUM_LIMIT:
        # GNU and BSD `sum` default to this algorithm and are much faster than
        # Python; a SysV `sum` reports 512-byte blocks and is not used.
        try:
            out = subprocess.run([exe, str(path)], capture_output=True, text=True, check=True).stdout.split()
            if len(out) >= 2 and out[0].isdigit() and out[1].isdigit() and int(out[1]) == blocks:
                return int(out[0]), blocks
        except (OSError, subprocess.CalledProcessError):
            pass
    return _bsd_sum_python(path)


def select_files(dna_checksums: T.Iterable[str], gtf_checksums: T.Iterable[str],
                 release: int) -> T.Tuple[str, str, str, str]:
    """Pick (sequence type, FASTA name, GTF name, assembly) from directory listings.

    The unmasked primary assembly is preferred; species without one (most
    non-human genomes) only have the toplevel FASTA."""
    dna_checksums, gtf_checksums = list(dna_checksums), list(gtf_checksums)
    for sequence in ("primary_assembly", "toplevel"):
        suffix = f".dna.{sequence}.fa.gz"
        names = sorted(n for n in dna_checksums if n.endswith(suffix))
        if names:
            break
    else:
        raise GenomeFetchError("No *.dna.primary_assembly.fa.gz or *.dna.toplevel.fa.gz file listed in CHECKSUMS.")
    fasta = names[0]
    prefix = fasta[:-len(suffix)]  # <Species>.<Assembly>
    gtf = f"{prefix}.{release}.gtf.gz"
    if gtf not in gtf_checksums:
        # The plain annotation is <Species>.<Assembly>.<N>.gtf.gz; variants
        # (.chr, .abinitio, .chr_patch_hapl_scaff) end differently. Worm, fly and
        # yeast in Ensembl use the Ensembl Genomes release as N (e.g. 63 in 116).
        pattern = re.compile(re.escape(prefix) + r"\.\d+\.gtf\.gz")
        names = [n for n in gtf_checksums if pattern.fullmatch(n)]
        if len(names) != 1:
            raise GenomeFetchError(f"Cannot identify the GTF for {prefix} among: {', '.join(gtf_checksums) or 'none'}.")
        gtf = names[0]
    assembly = prefix.split(".", 1)[1] if "." in prefix else ""
    return sequence, fasta, gtf, assembly


# ----------------------------------------------------------------- planning

def plan_fetch(species: str, release: T.Union[int, str, None] = "latest", division: str = "vertebrates",
               mirror: T.Optional[str] = None, rest_url: T.Optional[str] = ENSEMBL_REST) -> FetchPlan:
    """Resolve release, species and the exact FASTA/GTF URLs without downloading them."""
    release = resolve_release(release, division, mirror, rest_url)
    root = release_root(division, release, mirror)
    row = resolve_species(species, load_species_table(division, release, mirror), division, release, rest_url)
    path = _species_path(row)
    dna_dir, gtf_dir = f"{root}/fasta/{path}/dna", f"{root}/gtf/{path}"
    dna_text = _get_text(f"{dna_dir}/CHECKSUMS", missing_ok=True)
    if dna_text is None:
        raise GenomeFetchError(
            f"No genome FASTA for '{row['species']}' in Ensembl {division} release {release} "
            f"({dna_dir}/CHECKSUMS not found). Check the species (`uprobe genome search`), --division "
            "and --release; Ensembl Genomes divisions have their own release numbers."
        )
    gtf_text = _get_text(f"{gtf_dir}/CHECKSUMS", missing_ok=True)
    if gtf_text is None:
        raise GenomeFetchError(f"No GTF annotation for '{row['species']}' in Ensembl {division} "
                               f"release {release} ({gtf_dir}/CHECKSUMS not found).")
    dna_sums, gtf_sums = parse_checksums(dna_text), parse_checksums(gtf_text)
    sequence, fasta, gtf, assembly = select_files(dna_sums, gtf_sums, release)
    return FetchPlan(
        division=division, release=release, species=row["species"],
        display_name=row.get("name") or row["species"], assembly=assembly, sequence=sequence,
        fasta=RemoteFile("fasta", f"{dna_dir}/{fasta}", fasta, dna_sums[fasta]),
        gtf=RemoteFile("gtf", f"{gtf_dir}/{gtf}", gtf, gtf_sums[gtf]),
    )


def default_genome_name(plan: FetchPlan) -> str:
    """e.g. homo_sapiens_GRCh38_e116, arabidopsis_thaliana_TAIR10_eg63."""
    tag = "e" if plan.division == "vertebrates" else "eg"
    parts = [plan.species, plan.assembly, f"{tag}{plan.release}"]
    return re.sub(r"[^A-Za-z0-9._-]", "_", "_".join(p for p in parts if p))


# ---------------------------------------------------------------- download

def _progress(done: int, total: T.Optional[int]) -> str:
    if total:
        return f"{done / 2**20:.1f} / {total / 2**20:.1f} MiB ({100 * done / total:.0f}%)"
    return f"{done / 2**20:.1f} MiB"


def _download(url: str, dest: Path, size_hint: T.Optional[int] = None) -> None:
    """Download ``url`` to ``dest`` through ``dest.part``, resuming a partial file."""
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, _ATTEMPTS + 1):
        offset = part.stat().st_size if part.exists() else 0
        try:
            with _open(url, offset) as resp:
                if offset and resp.status != 206:
                    offset = 0  # the server ignored the Range header: start over
                if offset:
                    log.info(f"Resuming {dest.name} at {_progress(offset, None)}")
                length = resp.headers.get("Content-Length")
                total = offset + int(length) if length else size_hint
                done, last = offset, time.monotonic()
                with open(part, "ab" if offset else "wb") as fh:
                    while chunk := resp.read(_CHUNK):
                        fh.write(chunk)
                        done += len(chunk)
                        if time.monotonic() - last >= 10:
                            last = time.monotonic()
                            log.info(f"  {dest.name}: {_progress(done, total)}")
                log.info(f"  {dest.name}: {_progress(done, total)} downloaded")
            break
        except urllib.error.HTTPError as e:
            if e.code == 416 and offset:
                break  # the partial file is already complete; the checksum decides
            raise GenomeFetchError(f"Download failed: HTTP {e.code} for {url}") from e
        except (OSError, http.client.HTTPException) as e:
            if attempt == _ATTEMPTS:
                raise GenomeFetchError(
                    f"Download of {url} failed: {e}. Re-run the same command to resume from {part}."
                ) from e
            log.warning(f"Download interrupted ({e}); retrying ({attempt}/{_ATTEMPTS - 1})...")
            time.sleep(2 * attempt)
    os.replace(part, dest)


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def _fetch_file(remote: RemoteFile, directory: Path, record: T.Optional[dict],
                force: bool = False, keep_compressed: bool = False) -> dict:
    """Download, verify and decompress one file; return its manifest record."""
    target = directory / remote.local_name
    gz = directory / remote.name
    expected = _fmt_checksum(remote.checksum)
    if (not force and record and record.get("file") == remote.name and record.get("checksum") == expected
            and target.is_file() and target.stat().st_size == record.get("size")):
        log.info(f"{target.name}: already present and verified, skipping")
        return record
    if force:
        gz.unlink(missing_ok=True)
        gz.with_name(gz.name + ".part").unlink(missing_ok=True)
    for attempt in (1, 2):
        if gz.is_file():
            log.info(f"Found {gz.name}; verifying it instead of downloading")
        else:
            log.info(f"Downloading {remote.url} (~{remote.checksum[1] / 1024:.1f} MiB)")
            _download(remote.url, gz, remote.checksum[1] * 1024)
        actual = bsd_sum(gz)
        if actual == remote.checksum:
            log.info(f"{gz.name}: checksum OK ({expected}, {CHECKSUM_ALGORITHM})")
            break
        gz.unlink()
        msg = f"Checksum mismatch for {gz.name}: CHECKSUMS lists {expected}, the download has {_fmt_checksum(actual)}"
        if attempt == 2:
            raise GenomeFetchError(f"{msg}. The corrupt file was deleted; re-run the command to download it "
                                   "again, or try another --mirror.")
        log.warning(f"{msg}; downloading it again")
    # pyfaidx, makeblastdb and jellyfish need plain files (Ensembl's .gz is not BGZF).
    log.info(f"Decompressing {gz.name}")
    tmp = target.with_name(target.name + ".part")
    with gzip.open(gz, "rb") as src, open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst, _CHUNK)
    os.replace(tmp, target)
    if not keep_compressed:
        gz.unlink()
    return {"url": remote.url, "file": remote.name, "checksum": expected, "path": target.name,
            "size": target.stat().st_size, "fetched_at": _now()}


# ---------------------------------------------------------------- manifest

def _read_manifest(directory: Path) -> dict:
    path = directory / MANIFEST_NAME
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    except ValueError:
        log.warning(f"Ignoring unreadable {path}")
        return {}
    return data if isinstance(data, dict) else {}


def _write_atomic(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _needs_replace(directory: Path, manifest: dict, plan: FetchPlan, force: bool) -> bool:
    """True if ``directory`` holds another genome/release that ``--force`` replaces."""
    keys = ("division", "release", "species", "assembly", "sequence")
    if not manifest or all(manifest.get(k) == getattr(plan, k) for k in keys):
        return False
    held = (f"{manifest.get('species')} {manifest.get('assembly')} "
            f"(Ensembl {manifest.get('division')} release {manifest.get('release')}, {manifest.get('sequence')})")
    if not force:
        raise GenomeFetchError(f"{directory} already holds {held}. Choose another --name or --out, "
                               "or pass --force to replace it (its files and indices are deleted).")
    log.warning(f"--force: {held} in {directory} is replaced")
    return True


def _clear_directory(directory: Path, manifest: dict) -> None:
    """Delete a previous genome's files and the indices derived from them."""
    for record in (manifest.get("files") or {}).values():
        if isinstance(record, dict) and record.get("path"):
            (directory / record["path"]).unlink(missing_ok=True)
    for name in _DERIVED:
        path = directory / name
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink(missing_ok=True)
    for fai in directory.glob("*.fai"):
        fai.unlink()


# ------------------------------------------------------------ genomes.yaml

def load_genomes(path: T.Union[str, Path]) -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise GenomeFetchError(f"{path} must be a YAML mapping of genome names to entries.")
    return data


def update_genomes_yaml(path: T.Union[str, Path], name: str, entry: dict) -> bool:
    """Add or replace one entry, keeping the others. Return False if unchanged.

    A new entry is appended to the file text, so comments survive; replacing
    an existing entry rewrites the file."""
    path = Path(path)
    data = load_genomes(path)
    if data.get(name) == entry:
        return False
    text = path.read_text(encoding="utf-8") if path.exists() else ""
    new_text = None
    if name not in data and text.strip():
        block = yaml.safe_dump({name: entry}, sort_keys=False, default_flow_style=False, allow_unicode=True)
        candidate = text + ("" if text.endswith("\n") else "\n") + block
        try:
            if yaml.safe_load(candidate) == {**data, name: entry}:
                new_text = candidate
        except yaml.YAMLError:
            pass
    if new_text is None:
        data[name] = entry
        new_text = yaml.safe_dump(data, sort_keys=False, default_flow_style=False, allow_unicode=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_atomic(path, new_text)
    return True


def _make_entry(plan: FetchPlan, directory: Path, existing: T.Optional[dict],
                align_index: T.Optional[T.List[str]], jellyfish: T.Optional[bool]) -> dict:
    """genomes.yaml entry in the format UProbeAPI reads; user settings
    (align_index, jellyfish, extra keys) of an existing entry are kept."""
    entry = dict(existing) if isinstance(existing, dict) else {}
    division = "" if plan.division == "vertebrates" else f"{plan.division} "
    entry.update({
        "description": f"{plan.display_name} {plan.assembly} "
                       f"(Ensembl {division}release {plan.release}, {plan.sequence})",
        "species": plan.species,
        "fasta": str(directory / plan.fasta.local_name),
        "gtf": str(directory / plan.gtf.local_name),
        "out": str(directory),
    })
    if align_index is not None:
        entry["align_index"] = list(align_index)
    else:
        entry.setdefault("align_index", ["bowtie2"])
    if jellyfish is not None:
        entry["jellyfish"] = bool(jellyfish)
    else:
        # Default DNA-mode protocols compute kmer_count, which requires jellyfish.
        entry.setdefault("jellyfish", shutil.which("jellyfish") is not None)
    entry["source"] = {
        "provider": "Ensembl",
        "division": plan.division,
        "release": plan.release,
        "assembly": plan.assembly,
        "sequence": plan.sequence,
        "fasta_url": plan.fasta.url,
        "gtf_url": plan.gtf.url,
        "manifest": str(directory / MANIFEST_NAME),
    }
    return entry


# -------------------------------------------------------------- public API

def fetch_genome(species: str, release: T.Union[int, str, None] = "latest", division: str = "vertebrates",
                 out_dir: T.Union[str, Path] = "genomes", genomes_yaml: T.Union[str, Path] = "genomes.yaml",
                 name: T.Optional[str] = None, align_index: T.Optional[T.List[str]] = None,
                 jellyfish: T.Optional[bool] = None, build_index: bool = False, threads: int = 10,
                 keep_compressed: bool = False, force: bool = False, dry_run: bool = False,
                 mirror: T.Optional[str] = None, rest_url: T.Optional[str] = ENSEMBL_REST) -> dict:
    """Download an Ensembl genome FASTA + GTF into ``out_dir/<name>/``, verify
    them, and add/update ``name`` in ``genomes_yaml``.

    Re-running is cheap: files already present and verified (recorded in
    ``genome_source.json``) are skipped, and interrupted downloads resume.
    Returns ``{"name", "entry", "plan", "directory", "genomes_yaml", "changed"}``.
    """
    plan = plan_fetch(species, release, division, mirror, rest_url)
    name = name or default_genome_name(plan)
    if not _NAME_RE.fullmatch(name):
        raise GenomeFetchError(f"Invalid genome name '{name}': use letters, digits, '.', '_' or '-'.")
    directory = (Path(out_dir) / name).resolve()
    genomes_yaml = Path(genomes_yaml).resolve()
    manifest = _read_manifest(directory)
    replace = _needs_replace(directory, manifest, plan, force)
    existing = load_genomes(genomes_yaml).get(name)
    entry = _make_entry(plan, directory, existing, align_index, jellyfish)
    if isinstance(existing, dict) and existing.get("fasta") not in (None, entry["fasta"]):
        log.warning(f"Entry '{name}' in {genomes_yaml} points to {existing['fasta']}; "
                    "it is updated to the downloaded genome.")
    result = {"name": name, "entry": entry, "plan": plan, "directory": directory,
              "genomes_yaml": genomes_yaml, "changed": False}
    if dry_run:
        return result

    if replace:
        _clear_directory(directory, manifest)
        manifest = {}
    directory.mkdir(parents=True, exist_ok=True)
    need = sum((f.checksum[1] * 1024) * _SPACE_FACTOR[f.kind] for f in (plan.fasta, plan.gtf)
               if not (directory / f.local_name).is_file())
    free = shutil.disk_usage(directory).free
    if free < need:
        log.warning(f"Only {free / 2**30:.1f} GiB free in {directory}; roughly {need / 2**30:.1f} GiB "
                    "may be needed for the downloads and decompressed files (indices need more).")
    files = dict(manifest.get("files") or {})
    manifest = {
        "provider": "Ensembl", "division": plan.division, "release": plan.release,
        "species": plan.species, "display_name": plan.display_name, "assembly": plan.assembly,
        "sequence": plan.sequence, "checksum_algorithm": CHECKSUM_ALGORITHM, "files": files,
    }
    for remote in (plan.fasta, plan.gtf):
        files[remote.kind] = _fetch_file(remote, directory, files.get(remote.kind), force, keep_compressed)
        _write_atomic(directory / MANIFEST_NAME, json.dumps(manifest, indent=2) + "\n")
    result["changed"] = update_genomes_yaml(genomes_yaml, name, entry)
    log.info(f"Genome '{name}' is registered in {genomes_yaml}")

    if build_index:
        from .api import UProbeAPI
        try:
            api = UProbeAPI(protocol_config={"genome": name}, genomes_config=genomes_yaml, output_dir=directory)
            api.build_genome_index(threads=threads)
        except Exception as e:
            raise GenomeFetchError(
                f"Genome '{name}' was downloaded and registered, but building indices failed: {e}. "
                "Check that the aligners are installed, then re-run with --build-index or use `uprobe build-index`."
            ) from e
    return result


def search(query: str, division: str = "vertebrates", release: T.Union[int, str, None] = "latest",
           mirror: T.Optional[str] = None, rest_url: T.Optional[str] = ENSEMBL_REST,
           limit: int = 20) -> T.Tuple[int, T.List[dict]]:
    """Search the species list of a release; return (release, matching rows)."""
    release = resolve_release(release, division, mirror, rest_url)
    table = load_species_table(division, release, mirror)
    if not table:
        raise GenomeFetchError(f"No species list is published for Ensembl {division} release {release}.")
    return release, search_species(query, table, limit=limit)


def _index_status(fasta: T.Optional[Path], aligner: str) -> T.Dict[str, T.Optional[bool]]:
    from .tools import _blast_index_exists, _bowtie2_index_exists
    check = {"bowtie2": _bowtie2_index_exists, "blast": _blast_index_exists}.get(aligner)
    if fasta is None or check is None:
        return {"genome": None, "transcript": None}
    return {"genome": check(fasta.parent / f"{aligner}_genome" / fasta.stem),
            "transcript": check(fasta.parent / f"{aligner}_transcript" / fasta.stem)}


def describe_genomes(genomes_yaml: T.Union[str, Path]) -> T.List[dict]:
    """Summarise genomes.yaml entries, including whether files and indices exist."""
    rows = []
    for name, entry in load_genomes(genomes_yaml).items():
        entry = entry if isinstance(entry, dict) else {}
        source = entry.get("source") if isinstance(entry.get("source"), dict) else {}
        paths = {key: Path(entry[key]) if entry.get(key) else None for key in ("fasta", "gtf")}
        rows.append({
            "name": name,
            "species": entry.get("species"),
            "description": entry.get("description"),
            "fasta": str(paths["fasta"]) if paths["fasta"] else None,
            "fasta_exists": bool(paths["fasta"] and paths["fasta"].is_file()),
            "gtf": str(paths["gtf"]) if paths["gtf"] else None,
            "gtf_exists": bool(paths["gtf"] and paths["gtf"].is_file()),
            "indexes": {a: _index_status(paths["fasta"], a) for a in entry.get("align_index") or []},
            "jellyfish": bool(entry.get("jellyfish", False)),
            "source": source or None,
        })
    return rows
