"""``mapped_genes``, the transcriptome off-target search: settings and results.

A ``mapped_genes`` attribute aligns a sequence to the transcriptome index built
from ``transcript.fa`` (records ``{gene_id}_{transcript_id}``) and counts the
distinct genes it hits::

    target_mapped_genes:
      target: target_region
      type: mapped_genes
      aligner: bowtie       # bowtie2 (default heuristic) or bowtie (exhaustive)
      max_mismatches: 2     # bowtie only: every site within 0-3 substitutions
      strand: sense         # both (default), sense or antisense

``aligner: bowtie2`` runs Bowtie2 in ``--very-sensitive-local`` mode with
``-k 100``: local alignments may be partial, and as a seed-based heuristic it
can miss sites with mismatches. Its default minimum local score,
20 + 8 ln(length), exceeds the score of a perfect match (2 per base) below
22 nt, so shorter sequences never align and count 0 genes, and sequences
of a few nt more align only without mismatches. ``aligner: bowtie`` runs Bowtie 1 with
``-v <max_mismatches> -a``, which ignores base qualities and reports every
end-to-end alignment with at most that many substitutions (indels are not
considered), so the off-target sites within the limit are found exhaustively.

``strand`` selects the transcript strand a site may lie on, relative to the
aligned sequence: ``sense`` keeps sites where the sequence itself occurs in a
transcript (right for ``target_region``: its probes bind exactly these
windows), ``antisense`` sites where its reverse complement occurs (right for
a probe's binding part), and ``both`` either (the original behaviour).

Besides the count, the attribute lists the off-target genes of each sequence,
i.e. the genes other than the row's target, in ``<attribute>_off_targets``,
e.g. ``ACTG1:0, POTEF:2`` (gene name and fewest mismatches; Bowtie2 gives
names only). Genes are named after the GTF's gene_name; a name that several
gene ids share is written ``name(gene_id)``. With a ``gene_expression`` table
(see uprobe.core.gene_expression), off-target genes below its ``min_value``
(or below ``min_fraction_of_target`` of the target's value, or with
``absent: ignore`` missing from the table) are not counted and are listed in
``<attribute>_unexpressed`` instead. With ``exempt_copies`` (see
uprobe.core.copies), off-target genes that copy most of the target's sequence
are not counted either and are listed in ``<attribute>_copies``.
"""
import functools
import re
import shutil
from collections import defaultdict
from pathlib import Path

from .copies import copies_settings, find_copies
from .utils import get_logger

log = get_logger(__name__)

MAPPED_GENES_ALIGNERS = ('bowtie2', 'bowtie')
# Command-line options of both Bowtie and Bowtie2 for each strand setting.
STRAND_OPTIONS = {'both': [], 'sense': ['--norc'], 'antisense': ['--nofw']}
STRANDS = tuple(STRAND_OPTIONS)
DEFAULT_MAX_MISMATCHES = 2
MAX_MISMATCHES_LIMIT = 3  # bowtie -v accepts 0-3
# Shortest sequence Bowtie2's local mode (default --score-min G,20,8) can align.
BOWTIE2_LOCAL_MIN_LENGTH = 22
# Up to this length Bowtie2 misses many sites within 2 substitutions: on the
# Ensembl 113 human transcriptome it flagged 69-88% of 25-30-nt mRNA windows
# with such a site in another gene, and 97-99% of 40-52-nt windows.
BOWTIE2_INCOMPLETE_MAX_LENGTH = 30
# Columns a mapped_genes attribute adds next to its count.
OFF_TARGETS_SUFFIX = '_off_targets'
UNEXPRESSED_SUFFIX = '_unexpressed'
COPIES_SUFFIX = '_copies'
_MISSING_EXAMPLES = 20


def short_for_bowtie2(attr_name: str, short: int, total: int) -> str:
    """Message for sequences that a bowtie2 mapped_genes attribute cannot align."""
    return (f'{attr_name}: {short} of {total} sequences are shorter than {BOWTIE2_LOCAL_MIN_LENGTH} nt, which '
            f'Bowtie2\'s local mode cannot align (its minimum score exceeds a perfect match), so they count 0 '
            f'genes; use aligner: bowtie for an exhaustive search of short sequences')


def bowtie_unavailable(genome) -> str:
    """Why mapped_genes cannot search ``genome`` with aligner: bowtie; '' if it can.

    It can when the genome entry lists bowtie in align_index (the workflow
    then builds the Bowtie 1 transcriptome index if it is missing) and
    Bowtie 1 is installed: bowtie, and bowtie-build while there is no index.
    """
    if not isinstance(genome, dict):
        return 'the genome configuration is not known'
    if 'bowtie' not in (genome.get('align_index') or []):
        return "the genome's align_index has no bowtie"
    if not shutil.which('bowtie'):
        return 'bowtie (Bowtie 1) is not on PATH'
    if genome.get('fasta') and not shutil.which('bowtie-build'):
        from .tools import _bowtie_index_exists
        fasta = Path(genome['fasta'])
        if not _bowtie_index_exists(fasta.parent / 'bowtie_transcript' / fasta.stem):
            return 'its transcriptome index is not built and bowtie-build is not on PATH'
    return ''


def mapped_genes_settings(attr: dict, where: str = '') -> dict:
    """Return {aligner, strand, max_mismatches, exempt_copies} of a mapped_genes attribute.

    Raises ValueError starting with the option name (after ``where.``);
    ``max_mismatches`` is None for bowtie2, ``exempt_copies`` None when off.
    """
    prefix = f'{where}.' if where else ''
    aligner = attr.get('aligner')
    if aligner not in MAPPED_GENES_ALIGNERS:
        raise ValueError(f'{prefix}aligner must be one of {", ".join(MAPPED_GENES_ALIGNERS)}, got {aligner!r}')
    strand = attr.get('strand', 'both')
    if strand not in STRANDS:
        raise ValueError(f'{prefix}strand must be one of {", ".join(STRANDS)}, got {strand!r}')
    max_mismatches = None
    if aligner == 'bowtie':
        max_mismatches = attr.get('max_mismatches', DEFAULT_MAX_MISMATCHES)
        if type(max_mismatches) is not int or not 0 <= max_mismatches <= MAX_MISMATCHES_LIMIT:
            raise ValueError(f'{prefix}max_mismatches must be an integer from 0 to {MAX_MISMATCHES_LIMIT}, '
                             f'got {max_mismatches!r}')
    return {'aligner': aligner, 'strand': strand, 'max_mismatches': max_mismatches,
            'exempt_copies': copies_settings(attr.get('exempt_copies'), f'{prefix}exempt_copies')}


# -- gene names ------------------------------------------------------------

_GTF_ATTRIBUTES = {key: re.compile(rf'(?:^|;)\s*{key}\s+"([^"\r\n]*)"') for key in ('gene_id', 'gene_name', 'gene')}
_GTF_GENE_NUMBER = re.compile(r'GeneID:(\d+)')


def _gtf_attribute(info: str, key: str):
    match = _GTF_ATTRIBUTES[key].search(info)
    return (match.group(1).strip() or None) if match else None


@functools.lru_cache(maxsize=4)
def _read_gene_names(path: str, mtime: int, size: int) -> dict:
    """Gene key (as in transcript.fa records) -> gene name, from the GTF's exon lines.

    The gene id falls back as in gen.fun.read_gtf, which names the
    transcript.fa records; the name falls back to the `gene` attribute.
    """
    names = {}
    with open(path) as f:
        for line in f:
            fields = line.split('\t', 8)
            if len(fields) < 9 or fields[2] != 'exon' or line.startswith('#'):
                continue
            info = fields[8]
            gene_id = _gtf_attribute(info, 'gene_id')
            if gene_id is None:
                number = _GTF_GENE_NUMBER.search(info)
                gene_id = number.group(1) if number else (_gtf_attribute(info, 'gene_name')
                                                          or _gtf_attribute(info, 'gene'))
            if gene_id is None:
                continue
            key = gene_id.split('_')[0]  # as attributes.transcript_gene
            if not names.get(key):
                names[key] = _gtf_attribute(info, 'gene_name') or _gtf_attribute(info, 'gene')
    return names


class GeneAnnotation:
    """Names of the transcriptome's genes, to label off-targets and recognize the target."""

    def __init__(self, names: dict):
        self.names = names
        self.keys_by_name = defaultdict(set)
        for key, name in names.items():
            if name:
                self.keys_by_name[name].add(key)

    @classmethod
    def from_gtf(cls, gtf) -> 'GeneAnnotation':
        path = Path(gtf) if gtf else None
        if path is None or not path.is_file():
            return cls({})
        stat = path.stat()
        return cls(_read_gene_names(str(path.resolve()), stat.st_mtime_ns, stat.st_size))

    def label(self, key: str) -> str:
        name = self.names.get(key)
        if not name:
            return str(key)
        return name if len(self.keys_by_name[name]) == 1 else f'{name}({key})'

    def target_keys(self, target) -> set:
        """Genes of a target given by gene name (or gene id)."""
        return set(self.keys_by_name.get(str(target), ())) | {str(target)}


# -- off-target lists --------------------------------------------------------

def format_gene_hits(hits, annotation: GeneAnnotation) -> str:
    """'ACTG1:0, POTEF:2' from [(gene key, fewest mismatches or None)], fewest mismatches first."""
    labelled = sorted((-1 if mismatches is None else mismatches, annotation.label(key)) for key, mismatches in hits)
    return ', '.join(label if mismatches < 0 else f'{label}:{mismatches}' for mismatches, label in labelled)


def parse_gene_hits(text) -> list:
    """Gene labels of an <attribute>_off_targets cell, without mismatch counts."""
    if not isinstance(text, str):
        return []
    labels = []
    for item in text.split(','):
        item = item.strip()
        label, colon, mismatches = item.rpartition(':')
        if colon and mismatches.isdigit():
            item = label
        if item:
            labels.append(item)
    return labels


class OffTargetCounter:
    """Counts and lists the off-target genes of the mapped_genes attributes of one add_attributes call.

    The GTF's gene names and the gene_expression table are read once, when the
    first mapped_genes attribute needs them. ``report`` collects statistics
    for summary.json: {'gene_expression': {...}, 'attributes': {name: {...}}}.
    """

    def __init__(self, protocol: dict, genome: dict, report: dict = None):
        self.protocol = protocol
        self.genome = genome or {}
        self.report = report if report is not None else {}
        self.annotation = None
        self.expression = None

    def _load(self):
        if self.annotation is not None:
            return
        self.annotation = GeneAnnotation.from_gtf(self.genome.get('gtf'))
        if self.protocol.get('gene_expression') is None:
            return
        from .gene_expression import ExpressionTable
        table = ExpressionTable(self.protocol['gene_expression'])
        description = table.describe()
        names = self.annotation.names
        if names:
            found = sum(table.value(key, name) is not None for key, name in names.items())
            if found == 0:
                raise ValueError(f'gene_expression: none of the {len(names)} genes of the genome annotation is in '
                                 f'{description["table"]} (column {", ".join(table.gene_columns)}); check the gene '
                                 'identifiers and the species')
            description.update(annotation_genes=len(names), annotation_genes_in_table=found)
            relative = (f' or below {table.min_fraction_of_target:g} x the target\'s value'
                        if table.min_fraction_of_target is not None else '')
            log.info(f'gene_expression: {description["table"]} lists {found} of the {len(names)} annotated genes; '
                     f'off-target genes with {table.value_column} < {table.min_value:g}{relative} are not counted')
        self.report['gene_expression'] = description
        self.expression = table

    def _thresholds(self, target_keys: dict) -> dict:
        """Target -> the expression value an off-target gene of it must reach to count.

        That is min_value or, with min_fraction_of_target, that fraction of the
        target's value when it is higher. The target's value is the highest of
        its genes (a name can belong to several gene ids). Targets the table
        does not list keep min_value; summary.json lists them and the log warns.
        """
        expression = self.expression
        if expression.min_fraction_of_target is None:
            return {target: expression.min_value for target in target_keys}
        description = self.report['gene_expression']
        known = description.setdefault('targets', {})
        absent = description.setdefault('targets_not_in_table', [])
        thresholds, new_absent = {}, []
        for target, keys in target_keys.items():
            values = [expression.value(key, self.annotation.names.get(key)) for key in sorted(keys)]
            value = max((v for v in values if v is not None), default=None)
            thresholds[target] = expression.threshold(value)
            if not target:
                continue
            if value is not None:
                known[target] = {'value': value, 'threshold': thresholds[target]}
            elif target not in absent:
                absent.append(target)
                new_absent.append(target)
        absent.sort()
        if new_absent:
            log.warning(f'gene_expression: {len(new_absent)} target(s) not in the expression table, so only min_value '
                        f'applies to their off-target genes: {", ".join(sorted(new_absent)[:_MISSING_EXAMPLES])}')
        return thresholds

    def _transcript_fasta(self, aligner: str) -> Path:
        """transcript.fa of the <aligner>_transcript index, extracted again if it is gone."""
        from .tools import transcript_fasta
        fasta = Path(self.genome['fasta'])
        path = transcript_fasta(fasta, fasta.parent / f'{aligner}_transcript')
        if not path.is_file() or path.stat().st_size == 0:
            from .gen.fun import extract_trans_seqs
            log.info(f'exempt_copies needs the transcript sequences; extracting {path}')
            extract_trans_seqs(self.genome['gtf'], fasta, path)
        return path

    def _copies(self, df, pairs, target_keys: dict, sequence_hits: dict, settings: dict):
        """({(target, gene key): shared fraction} of the hit genes that copy their target, report) or ({}, None)."""
        options = settings.get('exempt_copies')
        if not options:
            return {}, None
        candidates = defaultdict(set)
        for sequence, target in pairs:
            candidates[target].update(key for key in sequence_hits.get(sequence, {}) if key not in target_keys[target])
        candidates = {target: genes for target, genes in candidates.items() if genes}
        if not candidates:
            return {}, {'targets_without_sequence': [], 'seconds': 0.0}
        return find_copies(df, candidates, self._transcript_fasta(settings['aligner']), options)

    def _copy_stats(self, attr_name: str, options: dict, copies: dict, copied: set, report: dict) -> dict:
        """summary.json entries of exempt_copies: the copies of each target with their shares; logs them."""
        label = self.annotation.label
        by_target = {}
        for (target, key), share in sorted(copies.items(), key=lambda item: (item[0][0], -item[1], label(item[0][1]))):
            by_target.setdefault(target, {})[label(key)] = round(share, 3)
        stats = {'exempt_copies': options, 'copy_genes': len(copied), 'copies': by_target,
                 'copies_seconds': report['seconds']}
        unknown = report['targets_without_sequence']
        if unknown:
            stats['targets_without_sequence'] = unknown
            log.warning(f'{attr_name}: exempt_copies cannot check {len(unknown)} target(s) whose windows give no '
                        f'sequence, so their off-target genes all count: {", ".join(unknown[:_MISSING_EXAMPLES])}')
        examples = '; '.join(f'{target}: ' + ', '.join(f'{name} ({share:.2f})' for name, share in list(genes.items())[:5])
                             for target, genes in list(by_target.items())[:5])
        log.info(f'{attr_name}: {len(copied)} off-target gene(s) contain at least {options["min_shared"]:g} of the '
                 f'{options["kmer"]}-mers of a target and are not counted for it (exempt_copies, '
                 f'{report["seconds"]:g} s)' + (f': {examples}' if examples else ''))
        return stats

    def annotate(self, df, attr_name: str, column: str, sequence_hits: dict, settings: dict) -> None:
        """Write <attr_name>, <attr_name>_off_targets and, as configured, _unexpressed and _copies.

        ``sequence_hits`` maps each sequence of ``df[column]`` to {gene key: fewest mismatches or None}.
        The count includes the target's own genes; another gene does not count
        when it copies the target (exempt_copies), when the expression table
        puts it below the target's threshold (see _thresholds) or, with
        gene_expression.absent: ignore, when the table does not list it.
        """
        self._load()
        annotation, expression = self.annotation, self.expression
        target_column = next((c for c in ('target', 'gene') if c in df.columns), None)
        targets = df[target_column].astype(str).tolist() if target_column else [''] * len(df)
        pairs = set(zip(df[column], targets))
        target_keys = {target: annotation.target_keys(target) for target in {target for _, target in pairs}}
        thresholds = self._thresholds(target_keys) if expression is not None else {}
        copies, copy_report = self._copies(df, pairs, target_keys, sequence_hits, settings)
        rows = {}
        hit, ignored, below_fraction, missing, copied = set(), set(), set(), set(), set()
        for sequence, target in pairs:
            count, counted, below, copy_hits = 0, [], [], []
            for key, mismatches in sequence_hits.get(sequence, {}).items():
                if key in target_keys[target]:
                    count += 1
                    continue
                hit.add(key)
                if (target, key) in copies:
                    copy_hits.append((key, mismatches))
                    copied.add(key)
                    continue
                value = expression.value(key, annotation.names.get(key)) if expression else None
                if expression is not None and value is None:
                    missing.add(key)
                    ruled_out = expression.absent == 'ignore'
                else:
                    ruled_out = value is not None and value < thresholds[target]
                    if ruled_out:
                        (ignored if value < expression.min_value else below_fraction).add(key)
                if ruled_out:
                    below.append((key, mismatches))
                    continue
                count += 1
                counted.append((key, mismatches))
            rows[(sequence, target)] = (count, format_gene_hits(counted, annotation),
                                        format_gene_hits(below, annotation), format_gene_hits(copy_hits, annotation))
        keys = list(zip(df[column], targets))
        df[attr_name] = [rows[key][0] for key in keys]
        df[attr_name + OFF_TARGETS_SUFFIX] = [rows[key][1] for key in keys]
        if expression is not None:
            df[attr_name + UNEXPRESSED_SUFFIX] = [rows[key][2] for key in keys]
        if copy_report is not None:
            df[attr_name + COPIES_SUFFIX] = [rows[key][3] for key in keys]
        stats = {'aligner': settings['aligner'], 'max_mismatches': settings['max_mismatches'],
                 'strand': settings['strand'], 'off_target_genes': len(hit)}
        if copy_report is not None:
            stats.update(self._copy_stats(attr_name, settings['exempt_copies'], copies, copied, copy_report))
        if expression is not None:
            labels = sorted(annotation.label(key) for key in missing)
            stats.update(below_min_value=len(ignored), missing_from_table=len(missing), absent=expression.absent,
                         missing_from_table_examples=labels[:_MISSING_EXAMPLES], missing_labels=labels)
            relative = ''
            if expression.min_fraction_of_target is not None:
                stats['below_target_fraction'] = len(below_fraction)
                relative = (f'; {len(below_fraction)} more are below gene_expression.min_fraction_of_target of a '
                            'target\'s value and not counted for it')
            absent = ('are not counted (gene_expression.absent: ignore)' if expression.absent == 'ignore'
                      else 'count as expressed')
            message = (f'{attr_name}: the candidates hit {len(hit)} off-target genes; {len(ignored)} are below '
                       f'gene_expression.min_value and not counted{relative}; {len(missing)} are not in the '
                       f'expression table and {absent}')
            if missing:
                log.warning(f'{message} (e.g. {", ".join(labels[:5])})')
            else:
                log.info(message)
        self.report.setdefault('attributes', {})[attr_name] = stats
