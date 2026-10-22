import sys
import logging
import os
import re
from os.path import exists
import typing as t
from pyfaidx import Fasta
from .sampling import iter_target_windows
import shutil
import subprocess

def get_logger(name):
    log = logging.getLogger(name)
    log.setLevel(logging.INFO)
    if not log.handlers:
        handler = logging.StreamHandler(sys.stderr)
        LOGGING_FMT = "%(name)-20s %(levelname)-7s @ %(asctime)s: %(message)s"
        LOGGING_DATE_FMT = "%m/%d/%y %H:%M:%S"
        handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
        log.addHandler(handler)
        log.setLevel(logging.DEBUG)
        log.propagate = False
    return log

def get_base_map():
    base_map = [b'\0' for i in range(256)]
    base_map[ ord('A') ] = b'T'
    base_map[ ord('T') ] = b'A'
    base_map[ ord('C') ] = b'G'
    base_map[ ord('G') ] = b'C'
    base_map[ ord('a') ] = b't'
    base_map[ ord('t') ] = b'a'
    base_map[ ord('c') ] = b'g'
    base_map[ ord('g') ] = b'c'
    base_map[ ord('N') ] = b'N'
    base_map[ ord('n') ] = b'n'
    base_map = bytes(b''.join(base_map))
    return base_map

BASEMAP = get_base_map()
def reverse_complement(seq):
    # Handle pandas Series by converting to string
    if hasattr(seq, 'iloc') and hasattr(seq, 'values'):  # Check if it's a pandas Series
        seq = str(seq.iloc[0]) if len(seq) > 0 else str(seq.values[0])
    elif hasattr(seq, '__iter__') and not isinstance(seq, str):
        # Handle other iterable types that might be passed
        seq = str(seq)
    
    res = seq[::-1]
    res = res.translate(BASEMAP)
    return res

def get_tmp_dir(basename):
    i = 0
    dirname = lambda: f"{basename}.{i}"
    while exists(dirname()):
        i += 1
    os.mkdir(dirname())
    return dirname()

def self_match(probe: str, min_match = 4):
    length = len(probe)
    probe_re = reverse_complement(probe)
    match_pairs = 0
    for i in range(0,length-min_match+1):
        tem = probe[i:min_match+i]
        for j in range(0,length-min_match+1):
            tem_re = probe_re[j:min_match+j]
            if tem == tem_re and i + j + min_match - 2 != length:
                match_pairs = match_pairs + 1
    return match_pairs

def write_fastq(outdir, gene, recname2seq: t.Mapping[str, str]):
    fq = f'{outdir}/{gene}.fq'
    with open(fq, 'w') as f:
        for recname, seq in recname2seq.items():
            f.write(f"@{recname}\n")
            f.write(seq+"\n")
            f.write("+\n")
            f.write("~"*len(seq)+"\n")
    return fq

# A genomic region target (source: genome) or avoid_otp target region.
REGION_PATTERN = re.compile(r'^(\S+):(\d+)-(\d+)$')


def parse_region(region: str) -> t.Tuple[str, int, int]:
    """(chrom, start, end) of 'chrom:start-end', 1-based and inclusive.

    As genome browsers and Ensembl write regions: chr1:1-100 is the first 100
    bases of chr1. Raises ValueError for anything else, a start below 1 or an
    end before the start.
    """
    match = REGION_PATTERN.match(str(region).strip())
    if not match:
        raise ValueError(f"region {region!r} is not chrom:start-end")
    chrom, start, end = match[1], int(match[2]), int(match[3])
    if start < 1:
        raise ValueError(f"region {region!r} starts at {start}; positions are 1-based, the first base is 1")
    if end < start:
        raise ValueError(f"region {region!r} ends before it starts")
    return chrom, start, end


# extracts.target_region.masked: what source: genome does with repeat-masked windows.
MASKED_OPTIONS = ('keep', 'drop')


def masked_prefix(seq: str):
    """Cumulative count of masked bases: soft-masked (lowercase) letters and N; element i counts seq[:i]."""
    import numpy as np
    codes = np.frombuffer(seq.encode('ascii', 'replace'), dtype=np.uint8)
    masked = ((codes >= ord('a')) & (codes <= ord('z'))) | (codes == ord('N'))
    return np.concatenate(([0], np.cumsum(masked)))


class RepeatMask:
    """Masked bases of the loci, from a file given apart from the genome (extracts.target_region.mask).

    A BED file (chrom, start, end: 0-based, end exclusive, in the coordinates of the genome FASTA;
    e.g. RepeatMasker or UCSC rmsk intervals, optionally gzipped), or a soft-masked FASTA of the
    loci whose record names are regions, chrom:start-end (1-based, inclusive, as targets are
    written), with repeats in lowercase. The intervals apply to the extracted locus, so the genome
    FASTA need not be a soft-masked copy.
    """

    def __init__(self, path):
        import numpy as np
        self.path = str(path)
        self.intervals = {}  # chrom -> (starts, ends) arrays, BED coordinates
        self.loci = {}       # chrom -> [(start, boolean array)] 1-based start of a soft-masked FASTA record
        opener = __import__('gzip').open if self.path.endswith('.gz') else open
        with opener(self.path, 'rt') as handle:
            first = next((line for line in handle if line.strip()), '')
        if first.startswith('>'):
            with Fasta(self.path) as fa:
                for name in fa.keys():
                    chrom, start, end = parse_region(name)
                    seq = fa[name][:].seq
                    if len(seq) != end - start + 1:
                        raise ValueError(f"{self.path}: record {name!r} has {len(seq)} bases, its name says "
                                         f"{end - start + 1}")
                    flags = np.frombuffer(seq.encode('ascii', 'replace'), dtype=np.uint8)
                    self.loci.setdefault(chrom, []).append(
                        (start, ((flags >= ord('a')) & (flags <= ord('z'))) | (flags == ord('N'))))
        else:
            rows = {}
            with opener(self.path, 'rt') as handle:
                for number, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith(('#', 'track', 'browser')):
                        continue
                    fields = line.split()
                    try:
                        start, end = int(fields[1]), int(fields[2])
                    except (IndexError, ValueError):
                        raise ValueError(f"{self.path}, line {number}: expected a BED line "
                                         f"'chrom start end', got {line.strip()[:60]!r}") from None
                    if start < 0 or end < start:
                        raise ValueError(f"{self.path}, line {number}: invalid interval {start}-{end}")
                    rows.setdefault(fields[0], []).append((start, end))
            for chrom, pairs in rows.items():
                array = np.array(pairs, dtype=np.int64)
                self.intervals[chrom] = (array[:, 0], array[:, 1])

    def chromosomes(self):
        return set(self.intervals) | set(self.loci)

    def flags(self, chrom, start, end):
        """Boolean array, one per base of chrom:start-end (1-based, inclusive): masked or not."""
        import numpy as np
        flags = np.zeros(end - start + 1, dtype=bool)
        if chrom in self.intervals:
            starts, ends = self.intervals[chrom]
            for a, b in zip(starts[(ends >= start) & (starts < end)], ends[(ends >= start) & (starts < end)]):
                flags[max(a + 1, start) - start:min(b, end) - start + 1] = True  # BED a..b is 1-based a+1..b
        for locus_start, locus in self.loci.get(chrom, []):
            low, high = max(locus_start, start), min(locus_start + len(locus) - 1, end)
            if low <= high:
                flags[low - start:high - start + 1] |= locus[low - locus_start:high - locus_start + 1]
        return flags


def extract_fasta(fasta_path, target, min_length, overlap=None, step=None, masked='keep', mask=None):
    """Windows of a region target, chrom:start-end (1-based, inclusive; see parse_region).

    sub_region numbers each window along the region, 1-based: window s-e
    covers the genomic bases start + s - 1 to start + e - 1. Windows are
    upper-cased. ``masked='drop'`` skips the windows with a masked base, i.e.
    a soft-masked (lowercase) base or an N of a repeat-masked FASTA, as
    OligoMiner mines a repeat-masked genome; the others keep the probe_id
    they have with ``masked='keep'``. ``mask`` (a RepeatMask) adds the bases it
    marks to the masked ones, for ``masked='drop'``, so an unmasked genome can be used.
    """
    if masked not in MASKED_OPTIONS:
        raise ValueError(f"masked must be one of {', '.join(MASKED_OPTIONS)}, got {masked!r}")
    chrom, start, end = parse_region(target)
    with Fasta(str(fasta_path)) as fa:
        raw = fa[chrom][start - 1:end].seq
    seq = raw.upper()
    prefix = masked_prefix(raw) if masked == 'drop' else None
    mask_marks = False  # the mask marks a base of the region
    if prefix is not None and mask is not None:
        import numpy as np
        flags = mask.flags(chrom, start, end)
        mask_marks = bool(flags.any())
        prefix = np.concatenate(([0], np.cumsum(((np.diff(prefix) > 0) | flags[:len(raw)]))))
    windows, skipped = [], 0
    for n, (i, stop, tem) in enumerate(iter_target_windows(seq, min_length, step, overlap), start=1):
        if prefix is not None and prefix[stop] != prefix[i]:
            skipped += 1
            continue
        windows.append([f"{target}_{n}", target, f"{i + 1}-{stop}", tem])
    if prefix is not None:
        log = get_logger(__name__)
        if raw == seq and mask is None:
            log.warning(f"{target} has no soft-masked (lowercase) base: if the FASTA is not repeat-masked, "
                        "masked: drop only skips windows with N (use a soft-masked genome, e.g. "
                        "`uprobe genome fetch --soft-masked`)")
        elif raw == seq and not mask_marks:
            elsewhere = (f"; the mask has no interval on {chrom} at all (its chromosomes: "
                         f"{', '.join(sorted(mask.chromosomes())[:5]) or 'none'}): check its chromosome names"
                         if chrom not in mask.chromosomes() else "")
            log.warning(f"{target} has no soft-masked (lowercase) base and the mask {mask.path} marks none of it"
                        f"{elsewhere}: masked: drop only skips windows with N here")
        log.info(f"{target}: masked: drop skipped {skipped} of {skipped + len(windows)} windows "
                 "with a soft-masked or N base")
    return windows

def gene_barcode(config: dict) -> dict:
    """
    Generates a dictionary of gene names to anchor barcodes
    """
    gene_barcode_dict = {}
    for target in config['targets']:
        if target in config['encoding']:
            barcodes = config['encoding'][target]
            # Get all barcodes for this target
            barcode_values = []
            for barcode_key in barcodes:
                barcode_value = config['barcode_set'][barcodes[barcode_key]]
                barcode_values.append(barcode_value)
            gene_barcode_dict[target] = tuple(barcode_values)
    return gene_barcode_dict

def check_and_install_tools(tools_list: t.List[str]) -> None:
    """
    Check if tools are available in PATH. If not, try to install via conda.
    """
    missing_tools = []
    for tool in tools_list:
        if not shutil.which(tool):
            missing_tools.append(tool)
    
    if not missing_tools:
        return

    log = get_logger(__name__)
    log.info(f"Missing tools: {', '.join(missing_tools)}. Attempting auto-installation via conda...")
    
    # Try installing all missing tools at once
    # Mapping tool names to package names if they differ
    # bowtie2 -> bowtie2
    # jellyfish -> jellyfish
    # blast -> blast
    # mmseqs -> mmseqs2
    package_map = {
        'mmseqs': 'mmseqs2'
    }
    
    packages = [package_map.get(tool, tool) for tool in missing_tools]
    
    cmd = ["conda", "install", "-y", "-c", "bioconda"] + packages
    
    try:
        # Check if conda is available
        if not shutil.which("conda"):
            raise EnvironmentError("Conda is not available in PATH. Cannot auto-install tools.")
            
        log.info(f"Running command: {' '.join(cmd)}")
        subprocess.check_call(cmd)
        log.info("Tools installed successfully.")
    except subprocess.CalledProcessError as e:
        log.error(f"Failed to install tools via conda: {e}")
        log.error("Please install the following tools manually: " + ", ".join(missing_tools))
        sys.exit(1)
    except Exception as e:
        log.error(f"Error during tool installation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    fasta_path = "/data/zhangqian/genomes/hg38/hg38.fa"
    target = "NC_000017.11:40304334-40362489"
    seq_list = extract_fasta(fasta_path, target)
    for id, region, seq in seq_list:
        print(id, region, seq)