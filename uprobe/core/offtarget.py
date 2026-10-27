"""``mapped_genes``, the transcriptome off-target search: settings and results.

A ``mapped_genes`` attribute aligns a sequence to the transcriptome index built
from ``transcript.fa`` (records ``{gene_id}_{transcript_id}``) and counts the
distinct genes it hits, the target counting once even when several gene ids
carry its name (copies on assembly patches or alternate haplotypes)::

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
``<attribute>_unexpressed`` instead. A table of transcripts (``level:
transcript``) values each off-target gene by the transcripts the sequence hits,
and ``<attribute>_hit_transcripts`` lists them with their values. With
``exempt_copies`` (see uprobe.core.copies), off-target genes that copy most of
the target's sequence are not counted either and are listed in
``<attribute>_copies``.
"""
import functools
import os
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from .copies import UNKNOWN_BIOTYPE, copies_settings, find_copies
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
HIT_TRANSCRIPTS_SUFFIX = '_hit_transcripts'
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
        from .genome_paths import index_prefix
        if not _bowtie_index_exists(index_prefix(genome, 'bowtie', 'transcript')):
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

_GTF_ATTRIBUTES = {key: re.compile(rf'(?:^|;)\s*{key}\s+"([^"\r\n]*)"')
                   for key in ('gene_id', 'gene_name', 'gene', 'transcript_id', 'gene_biotype', 'gene_type')}
_GTF_GENE_NUMBER = re.compile(r'GeneID:(\d+)')


def _gtf_attribute(info: str, key: str):
    match = _GTF_ATTRIBUTES[key].search(info)
    return (match.group(1).strip() or None) if match else None


@functools.lru_cache(maxsize=4)
def _read_genes(path: str, mtime: int, size: int, transcripts: bool = False) -> tuple:
    """(gene key -> gene name, gene key -> {transcript.fa record: transcript id}, gene key -> biotype)
    from the GTF's exon lines.

    Gene keys are as in transcript.fa records ``{gene_id}_{transcript_id}``:
    the gene id falls back as in gen.fun.read_gtf, which names the records,
    and the name falls back to the `gene` attribute. The biotype is the
    `gene_biotype` attribute (Ensembl) or `gene_type` (GENCODE); genes with
    neither are left out. The transcripts are read only when asked for
    (otherwise {}).
    """
    names, records, biotypes = {}, defaultdict(dict), {}
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
            if key not in biotypes:
                biotype = _gtf_attribute(info, 'gene_biotype') or _gtf_attribute(info, 'gene_type')
                if biotype:
                    biotypes[key] = biotype
            if transcripts:
                transcript_id = _gtf_attribute(info, 'transcript_id')
                if transcript_id is not None:
                    records[key][f'{gene_id}_{transcript_id}'] = transcript_id
    return names, dict(records), biotypes


class GeneAnnotation:
    """Names, biotypes (and, if read, transcripts) of the transcriptome's genes, to label off-targets,
    recognize the target and tell which copies exempt_copies may exempt."""

    def __init__(self, names: dict, transcripts: dict = None, biotypes: dict = None):
        self.names = names
        self.transcripts = transcripts or {}  # gene key -> {transcript.fa record: transcript id}
        self.biotypes = biotypes or {}  # gene key -> GTF gene_biotype / gene_type
        self.keys_by_name = defaultdict(set)
        for key, name in names.items():
            if name:
                self.keys_by_name[name].add(key)

    @classmethod
    def from_gtf(cls, gtf, transcripts: bool = False) -> 'GeneAnnotation':
        path = Path(gtf) if gtf else None
        if path is None or not path.is_file():
            return cls({})
        stat = path.stat()
        return cls(*_read_genes(str(path.resolve()), stat.st_mtime_ns, stat.st_size, transcripts))

    def transcript_id(self, record: str) -> str:
        """Transcript id of a transcript.fa record ``{gene_id}_{transcript_id}``."""
        found = self.transcripts.get(record.split('_')[0], {}).get(record)
        return found if found is not None else record.split('_', 1)[-1]

    def gene_transcripts(self, key: str) -> list:
        """Transcript ids of a gene."""
        return list(self.transcripts.get(key, {}).values())

    def label(self, key: str) -> str:
        name = self.names.get(key)
        if not name:
            return str(key)
        return name if len(self.keys_by_name[name]) == 1 else f'{name}({key})'

    def target_keys(self, target) -> set:
        """Genes of a target given by gene name (or gene id)."""
        return set(self.keys_by_name.get(str(target), ())) | {str(target)}


# -- an index that U-Probe did not build ---------------------------------------

# Tools that list the record names of an index (``-n``).
_INSPECT_TOOLS = {'bowtie': 'bowtie-inspect', 'bowtie2': 'bowtie2-inspect'}


def index_record_names(aligner: str, index_prefix) -> list:
    """Record names of a Bowtie or Bowtie2 index (``<aligner>-inspect -n``), each up to its first
    whitespace as the aligners report hits; None when the inspect tool is not on PATH."""
    tool = shutil.which(_INSPECT_TOOLS.get(aligner, ''))
    if not tool:
        return None
    files = sorted(Path(str(index_prefix)).parent.glob(Path(str(index_prefix)).name + '.*'))
    stamp = tuple((str(f), f.stat().st_mtime_ns, f.stat().st_size) for f in files)
    return list(_inspect_names(tool, str(index_prefix), stamp))


@functools.lru_cache(maxsize=8)
def _inspect_names(tool: str, prefix: str, stamp: tuple) -> tuple:
    result = subprocess.run([tool, '-n', prefix], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'{os.path.basename(tool)} -n {prefix} failed: {result.stderr.strip()[:300]}')
    return tuple(line.split()[0] for line in result.stdout.splitlines() if line.strip())


def record_gene_problem(names, genes, what: str, remedy: str) -> str:
    """Why the record names of a transcriptome index do not give the genes ``genes`` ('' if they do).

    mapped_genes reads a hit's gene from its record name, ``<gene_id>_<transcript_id>``, up to the
    first '_' (attributes.transcript_gene). A record named otherwise counts as a gene of its own (an
    isoform of the target becomes an off-target gene) or is merged with others (NM_000003.1 -> NM).
    """
    names = list(names)
    wrong = [name for name in names if name.split('_')[0] not in genes]
    if not wrong:
        return ''
    examples = ', '.join(f"{name} -> gene {name.split('_')[0]!r}" for name in wrong[:4])
    return (f'{what}: {len(wrong):,} of {len(names):,} records are not named <gene_id>_<transcript_id> after the '
            f"genome's annotation (e.g. {examples}). mapped_genes reads each hit's gene from the record name up to "
            f"the first '_', so it would count wrong off-target genes. {remedy}")


def transcript_index_problem(genome: dict, aligner: str, prefix) -> tuple:
    """(error, note) of the transcriptome index at ``prefix`` given under ``indexes``: error '' when its records
    name the genes of the genome's GTF; note why they were not checked (no GTF, no ``<aligner>-inspect``)."""
    from .genome_paths import index_root
    what = f'indexes.{aligner}_transcript ({prefix})'
    annotation = GeneAnnotation.from_gtf(genome.get('gtf'))
    if not annotation.names:
        return '', f'{what}: the GTF is not readable, so the record names of the index are not checked'
    names = index_record_names(aligner, prefix)
    if names is None:
        return '', (f'{what}: {aligner}-inspect is not on PATH, so the record names of the index are not checked; '
                    "they must be <gene_id>_<transcript_id> of the genome's GTF")
    remedy = (f"Remove indexes.{aligner}_transcript so that U-Probe builds the index from {index_root(genome)}/"
              "transcript.fa, or build yours from that file (its records are named so).")
    return record_gene_problem(names, set(annotation.names), what, remedy), ''


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


def format_hit_transcripts(genes, annotation: GeneAnnotation) -> str:
    """'ACTG1: ENST1=812.4; POTEF: ENST2=0, ENST3=NA' from [(gene key, mismatches, [(transcript, value)])].

    Genes come in the order of format_gene_hits; NA marks a transcript the table does not list.
    """
    labelled = sorted((-1 if mismatches is None else mismatches, annotation.label(key), values)
                      for key, mismatches, values in genes)
    return '; '.join(f'{label}: ' + ', '.join(f'{transcript}=' + ('NA' if value is None else f'{round(value, 3):g}')
                                               for transcript, value in sorted(values))
                     for _, label, values in labelled)


class OffTargetCounter:
    """Counts and lists the off-target genes of the mapped_genes attributes of one add_attributes call.

    The GTF's gene names (and, for a table of transcripts, the genes'
    transcripts) and the gene_expression table are read once, when the first
    mapped_genes attribute needs them. ``report`` collects statistics for
    summary.json: {'gene_expression': {...}, 'attributes': {name: {...}}}.
    """

    def __init__(self, protocol: dict, genome: dict, report: dict = None):
        self.protocol = protocol
        self.genome = genome or {}
        self.report = report if report is not None else {}
        self.annotation = None
        self.expression = None
        block = protocol.get('gene_expression')
        # A table of transcripts values each off-target gene by the transcripts a sequence hits.
        self.transcript_level = isinstance(block, dict) and block.get('level') == 'transcript'

    def _load(self):
        if self.annotation is not None:
            return
        self.annotation = GeneAnnotation.from_gtf(self.genome.get('gtf'), transcripts=self.transcript_level)
        if self.protocol.get('gene_expression') is None:
            return
        from .gene_expression import ExpressionTable
        table = ExpressionTable(self.protocol['gene_expression'])
        description = table.describe()
        names = self.annotation.names
        relative = (f' or below {table.min_fraction_of_target:g} x the target\'s value'
                    if table.min_fraction_of_target is not None else '')
        if table.level == 'transcript':
            found = self._check_transcripts(table, description)
            log.info(f'gene_expression: {description["table"]} lists {found} of the '
                     f'{description["annotation_transcripts"]} annotated transcripts; off-target genes whose hit '
                     f'transcripts sum to {table.value_column} < {table.min_value:g}{relative} are not counted')
        elif names:
            found = sum(table.value(key, name) is not None for key, name in names.items())
            if found == 0:
                hint = ('; its identifiers look like transcript ids: set gene_expression.level: transcript'
                        if table.identifiers == 'transcript' else '')
                raise ValueError(f'gene_expression: none of the {len(names)} genes of the genome annotation is in '
                                 f'{description["table"]} (column {", ".join(table.gene_columns)}); check the gene '
                                 f'identifiers and the species{hint}')
            description.update(annotation_genes=len(names), annotation_genes_in_table=found)
            log.info(f'gene_expression: {description["table"]} lists {found} of the {len(names)} annotated genes; '
                     f'off-target genes with {table.value_column} < {table.min_value:g}{relative} are not counted')
        self.report['gene_expression'] = description
        self.expression = table

    def _check_transcripts(self, table, description: dict) -> int:
        """How many annotated transcripts a table of transcripts lists; ValueError if none (or no annotation)."""
        transcripts = self.annotation.transcripts
        if not transcripts:
            gtf = self.genome.get('gtf')
            raise ValueError('gene_expression.level: transcript needs the genome\'s GTF, which assigns the '
                             'transcripts to genes' + (f'; {gtf} lists none' if gtf else ''))
        total = found = genes = 0
        for records in transcripts.values():
            listed = sum(table.transcript_value(transcript) is not None for transcript in records.values())
            total, found, genes = total + len(records), found + listed, genes + bool(listed)
        if found == 0:
            hint = ('; its identifiers look like gene ids: set gene_expression.level: gene'
                    if table.identifiers == 'gene' else '')
            raise ValueError(f'gene_expression: none of the {total} transcripts of the genome annotation is in '
                             f'{description["table"]} (column {table.transcript_column}); check the transcript '
                             f'identifiers and the species{hint}')
        description.update(annotation_transcripts=total, annotation_transcripts_in_table=found,
                           annotation_genes=len(transcripts), annotation_genes_in_table=genes)
        return found

    def _gene_value(self, key: str):
        """A gene's value: by id or name, or for a table of transcripts the sum over its listed transcripts."""
        if self.expression.level == 'transcript':
            return self.expression.transcripts_value(self.annotation.gene_transcripts(key))[0]
        return self.expression.value(key, self.annotation.names.get(key))

    def _thresholds(self, target_keys: dict) -> dict:
        """Target -> the expression value an off-target gene of it must reach to count.

        That is min_value or, with min_fraction_of_target, that fraction of the
        target's value when it is higher. The target's value is the highest of
        its genes (a name can belong to several gene ids); for a table of
        transcripts a gene's value is the sum over all its transcripts. Targets
        the table does not list keep min_value; summary.json lists them and the
        log warns.
        """
        expression = self.expression
        if expression.min_fraction_of_target is None:
            return {target: expression.min_value for target in target_keys}
        description = self.report['gene_expression']
        known = description.setdefault('targets', {})
        absent = description.setdefault('targets_not_in_table', [])
        thresholds, new_absent = {}, []
        for target, keys in target_keys.items():
            values = [self._gene_value(key) for key in sorted(keys)]
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
        from .genome_paths import index_root
        from .tools import transcript_fasta
        fasta = Path(self.genome['fasta'])
        root = index_root(self.genome)
        path = transcript_fasta(fasta, root / f'{aligner}_transcript', root)
        if not path.is_file() or path.stat().st_size == 0:
            from .gen.fun import extract_trans_seqs
            log.info(f'exempt_copies needs the transcript sequences; extracting {path}')
            extract_trans_seqs(self.genome['gtf'], fasta, path)
        return path

    def _copies(self, df, pairs, target_keys: dict, sequence_hits: dict, settings: dict):
        """({(target, gene key): shared fraction} of the hit genes that copy their target and are exempt, report)
        or ({}, None); see copies.find_copies."""
        options = settings.get('exempt_copies')
        if not options:
            return {}, None
        candidates = defaultdict(set)
        for sequence, target in pairs:
            candidates[target].update(key for key in sequence_hits.get(sequence, {}) if key not in target_keys[target])
        candidates = {target: genes for target, genes in candidates.items() if genes}
        if not candidates:
            return {}, {'not_exempt': {}, 'targets_without_sequence': [], 'seconds': 0.0}
        return find_copies(df, candidates, self._transcript_fasta(settings['aligner']), options,
                           self.annotation.biotypes)

    def _copy_stats(self, attr_name: str, options: dict, copies: dict, copied: set, report: dict) -> dict:
        """summary.json entries of exempt_copies: the exempt copies of each target with their shares, and the
        copies not exempt for their biotype (copies_not_exempt); logs them."""
        label = self.annotation.label
        by_target = {}
        for (target, key), share in sorted(copies.items(), key=lambda item: (item[0][0], -item[1], label(item[0][1]))):
            by_target.setdefault(target, {})[label(key)] = round(share, 3)
        stats = {'exempt_copies': options, 'copy_genes': len(copied), 'copies': by_target,
                 'copies_seconds': report['seconds']}
        not_exempt = {}
        for (target, key), (share, reason) in sorted(report['not_exempt'].items(),
                                                     key=lambda item: (item[0][0], -item[1][0], label(item[0][1]))):
            not_exempt.setdefault(target, {})[label(key)] = {'share': round(share, 3), 'not_exempt': reason}
        if not_exempt:
            stats['copies_not_exempt'] = not_exempt
            kept = '; '.join(f'{target}: ' + ', '.join(f'{name} ({entry["share"]:.2f}, {entry["not_exempt"]})'
                                                        for name, entry in list(genes.items())[:5])
                             for target, genes in list(not_exempt.items())[:5])
            log.info(f'{attr_name}: {len(report["not_exempt"])} copy pair(s) stay counted, as only copies that are '
                     f'not coding are exempt (exempt_copies): {kept}')
        untyped = {key for (_, key), (_, reason) in report['not_exempt'].items() if reason == UNKNOWN_BIOTYPE}
        if untyped:
            log.warning(f'{attr_name}: exempt_copies does not exempt {len(untyped)} copy gene(s) whose biotype the GTF '
                        f'does not give (no gene_biotype or gene_type), so they count: '
                        f'{", ".join(sorted(label(key) for key in untyped)[:_MISSING_EXAMPLES])}')
        unknown = report['targets_without_sequence']
        if unknown:
            stats['targets_without_sequence'] = unknown
            log.warning(f'{attr_name}: exempt_copies cannot check {len(unknown)} target(s) whose windows give no '
                        f'sequence, so their off-target genes all count: {", ".join(unknown[:_MISSING_EXAMPLES])}')
        examples = '; '.join(f'{target}: ' + ', '.join(f'{name} ({share:.2f})' for name, share in list(genes.items())[:5])
                             for target, genes in list(by_target.items())[:5])
        log.info(f'{attr_name}: {len(copied)} off-target gene(s) contain at least {options["min_shared"]:g} of the '
                 f'{options["kmer"]}-mers of a target, are not coding and are not counted for it '
                 f'(exempt_copies, {report["seconds"]:g} s)' + (f': {examples}' if examples else ''))
        return stats

    def _hit_values(self, records: dict) -> dict:
        """Gene key -> (summed value or None, [hit transcripts not in the table], [(transcript, value)]).

        ``records`` are the transcript.fa records ``{gene_id}_{transcript_id}``
        that one sequence hits; each gene is valued by the transcripts of it
        that the sequence hits.
        """
        by_gene = defaultdict(list)
        for record in records:
            by_gene[record.split('_')[0]].append(self.annotation.transcript_id(record))  # as transcript_gene
        values = {}
        for key, transcripts in by_gene.items():
            listed = [(transcript, self.expression.transcript_value(transcript)) for transcript in transcripts]
            known = [value for _, value in listed if value is not None]
            values[key] = (sum(known) if known else None,
                           [transcript for transcript, value in listed if value is None], listed)
        return values

    def annotate(self, df, attr_name: str, column: str, sequence_hits: dict, settings: dict,
                 transcript_hits: dict = None) -> None:
        """Write <attr_name>, <attr_name>_off_targets and, as configured, _unexpressed, _hit_transcripts and _copies.

        ``sequence_hits`` maps each sequence of ``df[column]`` to {gene key: fewest mismatches or None}.
        The count includes the target once, however many gene ids carry its
        name (e.g. copies of the gene on assembly patches or alternate
        haplotypes, which have gene ids of their own); another gene does not count
        when it copies the target and is not protein_coding (exempt_copies),
        when the expression table puts it below the target's threshold (see
        _thresholds) or, with
        gene_expression.absent: ignore, when the table does not list it.
        For a table of transcripts, ``transcript_hits`` maps each sequence to
        the transcript.fa records it hits, and an off-target gene's value is
        the sum over the hit transcripts of it; with absent: count, a hit
        transcript that the table does not list makes the gene count, with
        ignore it adds nothing.
        """
        self._load()
        annotation, expression = self.annotation, self.expression
        by_transcript = expression is not None and expression.level == 'transcript'
        if by_transcript and transcript_hits is None:
            raise ValueError(f'{attr_name}: gene_expression.level: transcript needs the transcripts each sequence hits')
        target_column = next((c for c in ('target', 'gene') if c in df.columns), None)
        targets = df[target_column].astype(str).tolist() if target_column else [''] * len(df)
        pairs = set(zip(df[column], targets))
        target_keys = {target: annotation.target_keys(target) for target in {target for _, target in pairs}}
        thresholds = self._thresholds(target_keys) if expression is not None else {}
        copies, copy_report = self._copies(df, pairs, target_keys, sequence_hits, settings)
        rows, hit_values = {}, {}
        hit, ignored, below_fraction, missing, copied = set(), set(), set(), set(), set()
        missing_transcripts = set()
        for sequence, target in pairs:
            count, counted, below, copy_hits, valued = 0, [], [], [], []
            if by_transcript and sequence not in hit_values:
                hit_values[sequence] = self._hit_values(transcript_hits.get(sequence, {}))
            own = sequence_hits.get(sequence, {}).keys() & target_keys[target]
            count += bool(own)  # the target is one gene, whatever number of gene ids carry its name
            for key, mismatches in sequence_hits.get(sequence, {}).items():
                if key in own:
                    continue
                hit.add(key)
                if (target, key) in copies:
                    copy_hits.append((key, mismatches))
                    copied.add(key)
                    continue
                if expression is None:
                    count += 1
                    counted.append((key, mismatches))
                    continue
                if by_transcript:
                    value, unlisted, listed = hit_values[sequence].get(key, (None, [], []))
                    valued.append((key, mismatches, listed))
                    missing_transcripts.update(unlisted)
                else:
                    value = expression.value(key, annotation.names.get(key))
                    unlisted = value is None
                if unlisted:
                    missing.add(key)
                if unlisted and expression.absent == 'count':
                    ruled_out = False
                elif value is None:
                    ruled_out = True  # absent: ignore, and the table lists none of it
                else:
                    ruled_out = value < thresholds[target]
                    if ruled_out:
                        (ignored if value < expression.min_value else below_fraction).add(key)
                if ruled_out:
                    below.append((key, mismatches))
                    continue
                count += 1
                counted.append((key, mismatches))
            rows[(sequence, target)] = (count, format_gene_hits(counted, annotation),
                                        format_gene_hits(below, annotation), format_gene_hits(copy_hits, annotation),
                                        format_hit_transcripts(valued, annotation))
        keys = list(zip(df[column], targets))
        df[attr_name] = [rows[key][0] for key in keys]
        df[attr_name + OFF_TARGETS_SUFFIX] = [rows[key][1] for key in keys]
        if expression is not None:
            df[attr_name + UNEXPRESSED_SUFFIX] = [rows[key][2] for key in keys]
        if by_transcript:
            df[attr_name + HIT_TRANSCRIPTS_SUFFIX] = [rows[key][4] for key in keys]
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
            if by_transcript:
                stats['missing_transcripts'] = len(missing_transcripts)
            relative = ''
            if expression.min_fraction_of_target is not None:
                stats['below_target_fraction'] = len(below_fraction)
                relative = (f'; {len(below_fraction)} more are below gene_expression.min_fraction_of_target of a '
                            'target\'s value and not counted for it')
            if by_transcript:
                absent = ('add nothing to their value (gene_expression.absent: ignore)' if expression.absent == 'ignore'
                          else 'make them count as expressed')
                missing_text = (f'{len(missing)} have hit transcripts that are not in the expression table '
                                f'({len(missing_transcripts)} transcripts), which {absent}')
                below_text = 'have hit transcripts summing below gene_expression.min_value and are not counted'
            else:
                absent = ('are not counted (gene_expression.absent: ignore)' if expression.absent == 'ignore'
                          else 'count as expressed')
                missing_text = f'{len(missing)} are not in the expression table and {absent}'
                below_text = 'are below gene_expression.min_value and not counted'
            message = (f'{attr_name}: the candidates hit {len(hit)} off-target genes; {len(ignored)} {below_text}'
                       f'{relative}; {missing_text}')
            if missing:
                log.warning(f'{message} (e.g. {", ".join(labels[:5])})')
            else:
                log.info(message)
        self.report.setdefault('attributes', {})[attr_name] = stats
