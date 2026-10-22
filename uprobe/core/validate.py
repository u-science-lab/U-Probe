"""Protocol checks behind ``uprobe validate``.

Every problem becomes a :class:`Finding` with a severity, a location path such as
``probes.probe_1.parts.arm5.expr``, a message and a concrete suggestion.
Errors are problems that make ``uprobe run`` fail, or make a configured step
silently do nothing (an unknown filter column, a misspelled post_process step);
warnings are suspicious but harmless.

The checks reuse the runtime validators (normalize_sampling, compile_layout,
validate_layout_references, check_expression, parse_condition, the probe DAG and
_validate_and_normalize_protocol) and turn their exceptions into findings.
Once the static checks pass, every probe is built for each target on synthetic
target regions with the safe expression evaluator, so expressions that fail at
run time are reported too.
"""
import ast
import copy
import difflib
import json
import logging
import re
import shutil
import string
import typing as T
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace

import yaml

from .gen.expression import CONTEXT_NAMES, PROBE_INDEX, check_expression, expression_references, probe_index_nodes
from .hybridization import (FIELDS as HYBRIDIZATION_FIELDS, TM_METHODS, default_tm_window, resolve_conditions,
                            tm_method, unused_fields,
                            validate_protocol_hybridization)
from .layout import compile_layout, validate_layout_references
from .offtarget import (BOWTIE2_INCOMPLETE_MAX_LENGTH, BOWTIE2_LOCAL_MIN_LENGTH, MAPPED_GENES_ALIGNERS,
                        MAX_MISMATCHES_LIMIT, STRANDS, mapped_genes_settings)
from .process import parse_condition
from .sampling import normalize_sampling
from .utils import MASKED_OPTIONS, REGION_PATTERN


# Sections required by _validate_and_normalize_protocol.
REQUIRED_KEYS = ('genome', 'targets', 'extracts', 'encoding', 'probes')
OPTIONAL_KEYS = (
    'name', 'description', 'attributes', 'post_process', 'summary',
    'barcode_generation',  # uprobe generate-barcodes -p
    'hybridization',  # hybridization conditions for Tm
    'gene_expression',  # expression table for off-target counting
    'meta',  # probe-template metadata
    'species', 'probe_type', 'probe_name', 'probetype',  # read by the HTTP task API
)
SOURCES = ('genome', 'mRNA', 'exon', 'CDS', 'UTR')
REGION_KEYS = ('source', 'length', 'step', 'overlap', 'layout', 'transcripts', 'masked')
PART_KEYS = ('template', 'parts', 'expr', 'length')
LEGACY_TARGET_NAMES = ('gene_id', 'gene', 'gene_name')
# Attribute types dispatched by attributes.add_attributes -> accepted aligners (the first is suggested).
ATTRIBUTE_TYPES = {
    'gc_content': None, 'annealing_temperature': None, 'length': None,
    'fold_score': None, 'self_match': None, 'blocks': None,
    'mapped_sites': ('bowtie2',), 'mapped_genes': MAPPED_GENES_ALIGNERS, 'kmer_count': ('jellyfish',),
    'difference': None,
}
ATTRIBUTE_ALIASES = {'n_mapped_genes': 'mapped_genes'}
ATTRIBUTE_KEYS = ('target', 'type', 'aligner', 'threads', 'kmer_len', 'size', 'hybridization',
                  'max_mismatches', 'strand', 'exempt_copies', 'tm_method', 'of')
POST_PROCESS_KEYS = ('filters', 'sorts', 'avoid_otp', 'remove_overlap', 'off_target_pool', 'equal_space', 'summary')
DNA_ONLY_STEPS = ('avoid_otp',)
SORT_KEYS = ('is_ascending', 'is_descending')
# Attribute types whose lower values are better; sorting them descending ranks the worst first.
LOWER_IS_BETTER = {'fold_score': 'most structured (fold_score is -MFE: higher = more structure)',
                   'self_match': 'most self-complementary', 'kmer_count': 'most repetitive',
                   'mapped_genes': 'least specific'}
OPERATORS = ('<=', '>=', '<', '>', '==')
OPERATOR_FIXES = {'=>': '>=', '=<': '<=', '=': '=='}
SUMMARY_KEYS = ('report_name', 'attributes', 'report')
SUMMARY_STATS = ('count', 'mean', 'std', 'min', 'max', 'median', 'q25', 'q75', 'range')
ALIGN_INDEXES = {'bowtie2': ('bowtie2-build', 'bowtie2'), 'bowtie': ('bowtie-build', 'bowtie'),
                 'blast': ('makeblastdb',), 'mmseqs': ('mmseqs',)}
# Bioconda package of each aligner.
ALIGN_PACKAGES = {'bowtie2': 'bowtie2', 'bowtie': 'bowtie   # Bowtie 1', 'blast': 'blast', 'mmseqs': 'mmseqs2'}
BOWTIE_INSTALL = 'conda install -c bioconda bowtie'
# Columns written by UProbeAPI.generate_target_seqs, and by avoid_otp.
TARGET_COLUMNS = ('probe_id', 'target', 'sub_region', 'start', 'end', 'target_start', 'target_end',
                  'exon_name', 'utr_name', 'transcript_names', 'n_trans', 'layout_id', 'target_region')
OTP_COLUMNS = ('on_target', 'off_target', 'target_ratio')
IUPAC = set('ACGTURYSWKMBDHVN')
REGION_RE = REGION_PATTERN
DRY_RUN_TARGETS = 200

SNIPPETS = {
    'genome': 'genome: hg38   # a top-level key of genomes.yaml',
    'targets': ('targets:\n'
                '  - CD4               # gene name (source: mRNA, exon, CDS or UTR)\n'
                '  - chr1:1000-2000    # region (source: genome)\n'
                '  - my_seq: ACGTACGT  # literal sequence'),
    'extracts': ('extracts:\n'
                 '  target_region:\n'
                 '    source: mRNA   # genome | mRNA | exon | CDS | UTR\n'
                 '    length: 40     # or an inclusive range: [38, 42]\n'
                 '    step: 20'),
    'encoding': 'encoding:\n  CD4:\n    BC1: ACGTACGTACGT',
    'probes': ('probes:\n'
               '  probe_1:\n'
               '    template: "{part1}{part2}"\n'
               '    parts:\n'
               '      part1:\n'
               '        expr: "rc(target_region[0:20])"\n'
               '      part2:\n'
               '        expr: "encoding[target][\'BC1\']"'),
    'attributes': 'attributes:\n  probe_1_tm:\n    target: probe_1\n    type: annealing_temperature',
    'post_process': 'post_process:\n  filters:\n    probe_1_tm:\n      condition: probe_1_tm >= 50 & probe_1_tm <= 65',
    'summary': 'summary:\n  report_name: rna_report\n  attributes: [probe_1_tm]',
}
LAYOUT_HINT = ('layout needs a template naming each part once and part lengths that can sum to '
               'target_region.length, e.g.\n'
               'layout:\n'
               '  template: "{left}{gap}{right}"\n'
               '  parts:\n'
               '    left: {length: [19, 25]}\n'
               '    gap: {length: [0, 2]}\n'
               '    right: {length: [19, 25]}')
EXPR_HINT = ('write one Python expression, e.g. rc(target_region[0:20]) or encoding[target][\'BC1\']; '
             'quote fixed sequences: expr: "\'ACGT\'"')
CONDITION_HINT = ('write "<attribute> <op> <number>" comparisons joined by &, e.g. "tm >= 50 & tm <= 60" '
                  '(ops: <=, >=, <, >, ==; plain decimal numbers)')
GENOME_SNIPPET = ('hg38:\n'
                  '  fasta: /path/to/hg38.fa\n'
                  '  gtf: /path/to/hg38.gtf\n'
                  '  align_index: [bowtie2]\n'
                  '  jellyfish: false')
HYBRIDIZATION_HINT = ('numbers >= 0 for any of ' + ', '.join(HYBRIDIZATION_FIELDS) + ', e.g.\n'
                      'hybridization:\n'
                      '  temperature_c: 37\n'
                      '  na_mM: 390          # 2x SSC\n'
                      '  formamide_pct: 30')
TRANSCRIPTS_HINT = ('transcripts: mane   # all (default): every transcript; canonical: Ensembl_canonical; mane: '
                    'MANE_Select, else Ensembl_canonical; or a list of transcript ids, e.g. [ENST00000646664]')
KMER_STRAND_HINT = ("strand: sense   # both (default): canonical k-mers, both strands; sense: the k-mers as "
                    "written on the FASTA's forward strand (OligoMiner, PaintSHOP); antisense: their reverse "
                    "complements")
MASKED_HINT = ('masked: drop   # skip source: genome windows with a soft-masked (lowercase) or N base; '
               'keep (default) upper-cases them')
COPIES_SNIPPET = ('exempt_copies:      # or exempt_copies: true for these defaults\n'
                  '  kmer: 30          # k-mer length, 12-32\n'
                  '  min_shared: 0.5   # share of the target\'s k-mers that one transcript of a copy contains')
OTP_SNIPPET = ('avoid_otp:\n'
               '  chr1:1000-2000:          # a target\n'
               '    target_regions: [chr1:1000-2000]\n'
               '    density_thresh: 0.00001')
GENE_EXPRESSION_SNIPPET = ('gene_expression:\n'
                           '  table: expression.tsv   # TSV/CSV with a gene id or name column and a numeric column\n'
                           '  level: gene             # default; transcript for a table of transcript ids\n'
                           '  gene_column: gene_id    # default: the first column (transcript_column with level: '
                           'transcript)\n'
                           '  value_column: TPM       # default: the only numeric column\n'
                           '  min_value: 1            # off-target genes below it are not counted\n'
                           '  min_fraction_of_target: 0.05   # optional: nor those below 5% of the target\'s value\n'
                           '  absent: count           # default; ignore: genes missing from the table do not count')


@dataclass
class Finding:
    severity: str  # 'error' or 'warning'
    path: str  # location inside the file, e.g. probes.probe_1.parts.part1.expr
    message: str
    suggestion: str
    file: T.Optional[str] = None
    line: T.Optional[int] = None
    parts: tuple = field(default=(), repr=False, compare=False)
    doc: str = field(default='protocol', repr=False, compare=False)

    def to_dict(self) -> dict:
        return {key: getattr(self, key) for key in ('severity', 'path', 'message', 'suggestion', 'file', 'line')}


def format_path(parts) -> str:
    text = ''
    for part in parts:
        if isinstance(part, int):
            text += f'[{part}]'
        elif re.fullmatch(r'[\w-]+', str(part)):
            text += f'.{part}' if text else str(part)
        else:
            text += f'[{str(part)!r}]'
    return text


def _closest(word, options) -> T.Optional[str]:
    options = [str(option) for option in options]
    word = str(word)
    for option in options:
        if option.lower() == word.lower():
            return option
    matches = difflib.get_close_matches(word, options, n=1, cutoff=0.6)
    return matches[0] if matches else None


def _did_you_mean(word, options, what='valid options') -> str:
    options = list(dict.fromkeys(str(option) for option in options))
    match = _closest(word, options)
    if match is not None:
        ratio = difflib.SequenceMatcher(None, str(word), match).ratio()
        ties = [m for m in difflib.get_close_matches(str(word), options, n=3, cutoff=0.6)
                if difflib.SequenceMatcher(None, str(word), m).ratio() == ratio]
        if len(ties) > 1 and match in ties:
            return 'did you mean one of ' + ', '.join(f"'{m}'" for m in sorted(ties)) + '?'
        return f"did you mean '{match}'?"
    if not options:
        return ''
    shown = ', '.join(options[:12]) + (', ...' if len(options) > 12 else '')
    return f'{what}: {shown}'


def _sampling_hint(message: str) -> T.Tuple[str, str]:
    if message.startswith('step'):
        return 'step', 'step: 20   # positive integer: distance between window starts'
    if 'overlap' in message:
        return 'overlap', 'replace overlap with step (distance between window starts), e.g. step: 20'
    return 'length', 'length: 40   # positive integer, or an inclusive range: [38, 42]'


def _failure_hint(message: str, lengths) -> str:
    """Suggestion for an ExpressionError raised while building a probe."""
    detail = message.split(': ', 1)[-1]
    if detail.startswith('index ') and 'out of range' in detail:
        return (f'an index is outside target_region, which is {lengths[0]}-{lengths[1]} nt long; '
                'keep indices below the minimum length or count from the end, e.g. target_region[-13:]')
    if detail.startswith('key '):
        return 'a lookup key is missing; check the encoding entries and part names'
    if detail.startswith(('cannot ', 'slice bounds')) or 'expects a sequence' in detail:
        return 'pass sequences (strings) to functions, e.g. rc(target_region[0:20])'
    return 'fix the expression so that it returns a nucleotide string'


def _test_sequence(size: int) -> str:
    return ('ACGTTGCAAGCTTCGATGCA' * (size // 20 + 1))[:size]


def _topological_order(nodes):
    """Return (dependencies-first order, cycles) for probe DAG nodes."""
    order, cycles, state = [], [], {}

    def visit(node, stack):
        state[node] = 1
        stack.append(node)
        for dep in node.deps:
            if state.get(dep) == 1:
                cycles.append(stack[stack.index(dep):])
            elif dep not in state:
                visit(dep, stack)
        stack.pop()
        state[node] = 2
        order.append(node)

    for node in nodes:
        if node not in state:
            visit(node, [])
    return order, cycles


def _fasta_lengths(fasta) -> T.Dict[str, int]:
    fai = Path(str(fasta) + '.fai')
    if fai.is_file():
        lengths = {}
        for line in fai.read_text().splitlines():
            fields = line.split('\t')
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
        return lengths
    from pyfaidx import Fasta
    with Fasta(str(fasta)) as fa:
        return {name: len(fa[name]) for name in fa.keys()}


class _Validator:
    def __init__(self, protocol, genomes=None, check_targets=False, files=None, lines=None):
        self.protocol = protocol
        self.genomes = genomes
        self.lookup_targets = check_targets
        self.files = files or {}
        self.lines = lines or {}
        self.findings: T.List[Finding] = []
        self.source = None
        self.layout = None  # compile_layout result
        self.lengths = None  # (min, max) length of extracted target regions
        self.genome_targets, self.direct_targets = [], {}
        self.target_names: T.List[str] = []
        self.target_paths: T.Dict[str, tuple] = {}
        self.nodes: T.Dict[str, T.Tuple[tuple, dict]] = {}  # probe/part name -> (path, config)
        self.exprs: T.List[T.Tuple[str, tuple, str]] = []
        self.unused_parts: T.List[T.Tuple[str, tuple]] = []
        self.lookups: T.List[tuple] = []  # encoding lookups: (path, target keys or None, barcode)
        self.referenced: T.Set[str] = set()
        self.dag_ok = False
        self.attributes: dict = {}
        self.coding_targets = None  # genome targets with CDS features in the GTF, once it was read
        self.built_lengths: T.Dict[str, int] = {}  # probe/part name -> shortest sequence in the dry run

    # -- bookkeeping -----------------------------------------------------
    def add(self, severity, path, message, suggestion='', doc='protocol'):
        path = tuple(path)
        lines = self.lines.get(doc, {})
        line = next((lines[path[:i]] for i in range(len(path), 0, -1) if path[:i] in lines), None)
        finding = Finding(severity, format_path(path), message, suggestion,
                          self.files.get(doc, doc), line, path, doc)
        if finding not in self.findings:
            self.findings.append(finding)

    def error(self, path, message, suggestion='', doc='protocol'):
        self.add('error', path, message, suggestion, doc)

    def warning(self, path, message, suggestion='', doc='protocol'):
        self.add('warning', path, message, suggestion, doc)

    def has_errors(self, *sections) -> bool:
        return any(f.severity == 'error' and f.doc == 'protocol' and f.parts[:1] and f.parts[0] in sections
                   for f in self.findings)

    def region(self) -> dict:
        extracts = self.protocol.get('extracts')
        region = extracts.get('target_region') if isinstance(extracts, dict) else None
        return region if isinstance(region, dict) else {}

    def layout_parts(self) -> T.List[str]:
        layout = self.region().get('layout')
        return list(layout['parts']) if isinstance(layout, dict) and isinstance(layout.get('parts'), dict) else []

    def sequence_columns(self) -> T.List[str]:
        return ['target_region'] + list(self.nodes) + [f'target_parts.{name}' for name in self.layout_parts()]

    def known_columns(self, protocol) -> T.List[str]:
        """Columns available to filters, sorts and summary; attributes first."""
        columns = list(self.attributes)
        for name, attr in self.attributes.items():
            if isinstance(attr, dict) and str(attr.get('type', '')).lower() == 'mapped_sites':
                columns.append(f'{name}_num')
        columns += list(TARGET_COLUMNS) + self.sequence_columns()
        columns += [f'target_parts.{name}.{key}' for name in self.layout_parts() for key in ('length', 'start', 'end')]
        post_process = protocol.get('post_process')
        if isinstance(post_process, dict) and post_process.get('avoid_otp'):
            columns += list(OTP_COLUMNS)
        return list(dict.fromkeys(columns))

    # -- driver ----------------------------------------------------------
    def run(self) -> T.List[Finding]:
        if not isinstance(self.protocol, dict):
            self.error((), f'the protocol must be a YAML mapping, got {type(self.protocol).__name__}',
                       'start from the sections: ' + ', '.join(REQUIRED_KEYS))
            return self.findings
        self.check_top_level()
        self.check_extracts()
        self.check_targets()
        self.check_probes()
        self.check_encoding()
        self.check_hybridization()
        protocol, generated = self.normalize()
        self.check_attributes(protocol, 'attributes' in generated)
        self.check_gene_expression(protocol)
        if 'post_process' in generated:
            self.check_default_tm_filters(protocol)
        self.check_post_process(protocol)
        self.check_summary(protocol)
        if self.genomes is not None:
            self.check_genome(protocol, 'attributes' in generated)
        self.check_exon_source()
        self.dry_run()
        self.check_alignment_lengths(protocol, 'attributes' in generated)
        return self.findings

    # -- top level -------------------------------------------------------
    def check_top_level(self):
        p = self.protocol
        for key in p:
            if key not in REQUIRED_KEYS + OPTIONAL_KEYS:
                self.warning((key,), f'unknown top-level key {key!r} is ignored',
                             _did_you_mean(key, REQUIRED_KEYS + OPTIONAL_KEYS, 'known sections'))
        for key in REQUIRED_KEYS:
            if key not in p:
                self.error((key,), f'missing required section `{key}`', 'add it, e.g.\n' + SNIPPETS[key])
            elif key != 'encoding' and p[key] in (None, '', [], {}):
                self.error((key,), f'`{key}` is empty', 'fill it in, e.g.\n' + SNIPPETS[key])
        kinds = {'targets': (list, dict), 'extracts': dict, 'encoding': dict, 'probes': dict,
                 'attributes': dict, 'post_process': dict, 'summary': dict}
        for key, kind in kinds.items():
            value = p.get(key)
            if value not in (None, '') and not isinstance(value, kind):
                expected = 'a list' if kind == (list, dict) else 'a mapping'
                extra = ' (list-style probes are not supported; key probes by name)' if key == 'probes' else ''
                self.error((key,), f'`{key}` must be {expected}, got {type(value).__name__}{extra}',
                           'e.g.\n' + SNIPPETS[key])
        if p.get('genome') not in (None, '', [], {}) and not isinstance(p['genome'], str):
            self.error(('genome',), f'genome must be the name of a genome, got {p["genome"]!r}', SNIPPETS['genome'])
        if 'name' in p and (not isinstance(p['name'], str) or not p['name'].strip()):
            self.warning(('name',), f'name is {p["name"]!r}; output files are named after it', 'name: my_probe_design')

    # -- extracts --------------------------------------------------------
    def check_extracts(self):
        extracts = self.protocol.get('extracts')
        if not isinstance(extracts, dict) or not extracts:
            return
        for key in extracts:
            if key != 'target_region':
                self.warning(('extracts', key), f'unknown key {key!r} in extracts is ignored',
                             'only extracts.target_region is used')
        path = ('extracts', 'target_region')
        region = extracts.get('target_region')
        if not isinstance(region, dict):
            message = 'missing `extracts.target_region`' if region is None else 'target_region must be a mapping'
            self.error(path, message, 'e.g.\n' + SNIPPETS['extracts'])
            return
        for key in region:
            if key not in REGION_KEYS:
                self.warning(path + (key,), f'unknown key {key!r} in target_region is ignored',
                             _did_you_mean(key, REGION_KEYS, 'known keys'))
        source = region.get('source')
        if source is None:
            self.error(path + ('source',), 'missing `source`',
                       'source: one of genome (targets are regions such as chr1:100-200), '
                       'mRNA, exon, CDS or UTR (targets are gene names)')
        elif source not in SOURCES:
            self.error(path + ('source',), f'unknown source {source!r}', _did_you_mean(source, SOURCES))
        else:
            self.source = source
        self.check_transcripts(region, path)
        self.check_masked(region, path)
        if region.get('length') is None:
            self.error(path + ('length',), 'missing `length`',
                       'length: 40   # positive integer, or an inclusive range: [38, 42]')
            return
        try:
            lower, upper, _ = normalize_sampling(region['length'], region.get('step'), region.get('overlap'))
        except ValueError as exc:
            key, hint = _sampling_hint(str(exc))
            self.error(path + (key,), str(exc), hint)
            return
        if region.get('step') is not None and region.get('overlap') is not None:
            self.warning(path + ('overlap',), 'overlap is ignored because step is set',
                         'remove overlap (step is the distance between window starts)')
        self.lengths = (lower, upper)
        if 'layout' in region:
            try:
                self.layout = compile_layout(region['layout'], region['length'])
            except Exception as exc:
                self.error(path + ('layout',), str(exc), LAYOUT_HINT)
                return
            self.lengths = (min(self.layout[1]), max(self.layout[1]))

    def check_transcripts(self, region, path):
        from .gen.fun import TRANSCRIPT_RULES, transcripts_setting
        if 'transcripts' not in region:
            return
        value = region['transcripts']
        try:
            setting = transcripts_setting(value)
        except ValueError as exc:
            hint = _did_you_mean(value, TRANSCRIPT_RULES) if isinstance(value, str) else ''
            hint = hint + '\n' if hint.startswith('did') else ''
            self.error(path + ('transcripts',), str(exc), hint + TRANSCRIPTS_HINT)
            return
        if self.source == 'genome' and setting != 'all':
            self.warning(path + ('transcripts',), 'transcripts applies to gene targets (source mRNA, exon, CDS or '
                         'UTR) and is ignored for source: genome', 'remove transcripts')

    def check_masked(self, region, path):
        if 'masked' not in region:
            return
        value = region['masked']
        if value not in MASKED_OPTIONS:
            self.error(path + ('masked',), f'masked must be one of {", ".join(MASKED_OPTIONS)}, got {value!r}',
                       MASKED_HINT)
        elif value != 'keep' and self.source is not None and self.source != 'genome':
            self.warning(path + ('masked',), 'masked applies to source: genome and is ignored for gene targets',
                         'remove masked')

    # -- targets ---------------------------------------------------------
    def check_targets(self):
        targets = self.protocol.get('targets')
        if not isinstance(targets, (list, dict)) or not targets:
            return
        items = []
        if isinstance(targets, dict):
            items = [(('targets', name), name, seq) for name, seq in targets.items()]
        else:
            for i, item in enumerate(targets):
                if isinstance(item, str):
                    items.append((('targets', i), item, None))
                elif isinstance(item, dict):
                    items.extend((('targets', i, name), name, seq) for name, seq in item.items())
                else:
                    self.error(('targets', i), f'target entry {item!r} ({type(item).__name__}) is ignored',
                               'write a gene name or region as a string (quote values YAML reads as numbers), '
                               'or `name: SEQUENCE`')
        from .api import UProbeAPI
        self.genome_targets, self.direct_targets = UProbeAPI._parse_targets(SimpleNamespace(protocol=self.protocol))
        self.target_names = list(dict.fromkeys(list(self.genome_targets) + list(self.direct_targets)))
        seen = set()
        for path, name, seq in items:
            if name in seen:
                self.warning(path, f'duplicate target {name!r}', 'remove the duplicate entry')
            seen.add(name)
            self.target_paths.setdefault(name, path)
            if seq is not None and not isinstance(seq, str):
                self.error(path, f'the sequence of target {name!r} must be a string',
                           'write `name: ACGT...`, or list the name alone for a genome target')
            elif seq:
                bad = ''.join(sorted(set(seq.upper()) - IUPAC))
                if bad:
                    self.warning(path, f'the sequence of target {name!r} contains non-nucleotide characters {bad!r}',
                                 'use only A, C, G, T/U and IUPAC codes such as N')
                if self.lengths and len(seq) < self.lengths[0]:
                    self.error(path, f'the sequence of target {name!r} is {len(seq)} nt, shorter than the minimum '
                               f'target_region length {self.lengths[0]}, so it yields no candidates',
                               f'provide at least {self.lengths[0]} nt, or lower extracts.target_region.length')
            elif self.source == 'genome':
                match = REGION_RE.match(str(name))
                if not match:
                    self.error(path, f'target {name!r} is not a region; source: genome needs chrom:start-end',
                               'e.g. chr1:100000-101000, or use source: mRNA, exon, CDS or UTR for gene names')
                elif int(match[2]) < 1:
                    self.error(path, f'region {name!r} starts at {match[2]}; regions are 1-based and inclusive, '
                               'the first base is 1', f'e.g. {match[1]}:1-{match[3]}')
                elif int(match[2]) > int(match[3]):
                    self.error(path, f'region {name!r} ends before it starts', 'write chrom:start-end with start <= end')
            elif self.source is not None and REGION_RE.match(str(name)):
                self.warning(path, f'target {name!r} looks like a genomic region, but source is {self.source} '
                             '(gene names)', 'use extracts.target_region.source: genome for regions')

    # -- probes ----------------------------------------------------------
    def check_probes(self):
        probes = self.protocol.get('probes')
        if not isinstance(probes, dict) or not probes:
            return
        if 'layout' in self.region() and 'target_parts' in probes:
            self.error(('probes', 'target_parts'), 'the probe name target_parts is reserved when a layout is set',
                       'rename this probe')
        for name, config in probes.items():
            if name in ('target_region', 'target', 'encoding', 'rc'):
                self.warning(('probes', name), f'probe name {name!r} hides the expression name of the same name',
                             f'rename it, e.g. {name}_probe')
            self._walk_probe(str(name), name, config, ('probes', name), top=True)
        self._check_expressions()
        for full, path in self.unused_parts:
            if full not in self.referenced:
                part = full.rsplit('.', 1)[-1]
                self.warning(path, f'part {part!r} is not used in the template',
                             f'add {{{part}}} to the template, or remove the part')
        if not self.has_errors('probes'):
            self._check_cycles()

    def _walk_probe(self, full, key, config, path, top):
        kind = 'probe' if top else 'part'
        if re.search(r'[.:]', str(key)):
            self.error(path, f'{kind} name {key!r} must not contain "." or ":" (they separate probe and part names)',
                       f"rename it, e.g. {re.sub(r'[.:]', '_', str(key))}")
        if not isinstance(config, dict):
            self.error(path, f'{kind} {full!r} must be a mapping',
                       'e.g.\n' + (SNIPPETS['probes'] if top else 'expr: "rc(target_region[0:20])"'))
            return
        for option in config:
            if option not in PART_KEYS:
                hint = ('for a fixed sequence write expr: "\'ACGT\'"' if option in ('sequence', 'seq')
                        else _did_you_mean(option, PART_KEYS, 'known keys'))
                self.warning(path + (option,), f'unknown key {option!r} is ignored', hint)
        has_template, has_expr = 'template' in config, 'expr' in config
        if top and not has_template:
            if has_expr:
                self.error(path, f'top-level probe {full!r} needs `template` and `parts`; `expr` is only allowed '
                           'inside parts', 'wrap it:\ntemplate: "{part1}"\nparts:\n  part1:\n    expr: '
                           + json.dumps(str(config['expr'])))
            else:
                self.error(path, f'probe {full!r} needs `template` and `parts`', 'e.g.\n' + SNIPPETS['probes'])
            return
        if not has_template and not has_expr:
            self.error(path, f'part {full!r} needs `expr`, or `template` with `parts`',
                       'e.g. expr: "rc(target_region[0:20])"')
            return
        self.nodes[full] = (path + ('template' if has_template else 'expr',), config)
        if not has_template:
            if 'parts' in config:
                self.warning(path + ('parts',), '`parts` is ignored without `template`',
                             'add a template that uses the parts, or remove them')
            expr = config['expr']
            if not isinstance(expr, str) or not expr.strip():
                self.error(path + ('expr',), 'expr must be a non-empty string', EXPR_HINT)
            else:
                self.exprs.append((full, path + ('expr',), expr))
            return
        if has_expr:
            self.warning(path + ('expr',), '`expr` is ignored because `template` is set',
                         'keep either template + parts or expr')
        parts = config.get('parts')
        if not isinstance(parts, dict):
            self.error(path + ('parts',), f'`parts` must be a mapping of part name -> {{expr: ...}}, '
                       f'got {type(parts).__name__}', 'e.g.\nparts:\n  part1:\n    expr: "rc(target_region[0:20])"')
            parts = {}
        template = config['template']
        if not isinstance(template, str):
            self.error(path + ('template',), 'template must be a string', 'quote it: template: "{part1}{part2}"')
        else:
            self._check_template(full, path, template, parts)
        for part_name, part_config in parts.items():
            self._walk_probe(f'{full}.{part_name}', part_name, part_config, path + ('parts', part_name), top=False)

    def _check_template(self, full, path, template, parts):
        tpath = path + ('template',)
        try:
            fields = list(string.Formatter().parse(template))
        except ValueError as exc:
            self.error(tpath, f'invalid template: {exc}', 'use {part_name} placeholders; write literal braces as {{ }}')
            return
        names = [str(name) for name in parts]
        used = set()
        for literal, name, spec, conversion in fields:
            bad = ''.join(sorted(set(literal.upper()) - IUPAC))
            if bad:
                self.warning(tpath, f'literal text {literal!r} contains non-nucleotide characters {bad!r}',
                             'templates hold {placeholders} and bases only (A, C, G, T, N ...)')
            if name is None:
                continue
            if name == '' or name.isdigit():
                self.error(tpath, f'positional placeholder {{{name}}} is not supported',
                           'name each placeholder after a part, e.g. {part1}')
            elif spec or conversion:
                self.error(tpath, f'placeholder {{{name}}} uses a format spec or conversion',
                           f'use a plain placeholder: {{{name}}}')
            elif name not in names:
                hint = _did_you_mean(name, names, 'parts') if names else f'add parts.{name} with an expr'
                self.error(tpath, f'placeholder {{{name}}} has no matching part', hint)
            used.add(name)
        for name in names:
            if name not in used:
                self.unused_parts.append((f'{full}.{name}', path + ('parts', name)))

    def _check_expressions(self):
        names = set(self.nodes)
        known = names | set(CONTEXT_NAMES)
        region, layout_parts = self.region(), self.layout_parts()
        for node, path, expr in self.exprs:
            # The grammar of the safe evaluator that builds the probes; names
            # are checked below so that unknown ones get a suggestion.
            problems = check_expression(expr)
            for problem in problems:
                self.error(path, problem, EXPR_HINT)
            if problems:
                continue
            for ref in sorted(expression_references(expr)):
                head = ref.split('.')[0]
                if ref in names:
                    self.referenced.add(ref)
                elif ref in known:
                    continue
                elif head in names:
                    options = [n for n in names if n == head or n.startswith(head + '.')]
                    self.error(path, f'{ref!r} is not a probe or part', _did_you_mean(ref, options, f'parts of {head}'))
                elif '.' in ref:
                    self.error(path, f'unknown name {ref!r}: dotted names must refer to a probe part '
                               '(attribute access is not allowed)', EXPR_HINT)
                elif head in LEGACY_TARGET_NAMES:
                    fixed = re.sub(rf'\b{re.escape(head)}\b', 'target', expr)
                    self.error(path, f'unknown name {head!r}; the current target name is `target`',
                               f'replace {head} with target: expr: {json.dumps(fixed)}')
                elif re.fullmatch(r'[ACGTUNacgtun]+', head):
                    self.error(path, f'unknown name {head!r}', f'quote literal sequences: expr: "\'{head}\'"')
                else:
                    self.error(path, f'unknown name {head!r}', _did_you_mean(head, sorted(known), 'available names'))
            try:
                validate_layout_references({'extracts': {'target_region': region}, 'probes': {node: {'expr': expr}}})
            except Exception as exc:
                missing = re.search(r'missing target part (.+)$', str(exc))
                if not layout_parts:
                    hint = 'define extracts.target_region.layout with named parts, or slice target_region instead'
                elif missing:
                    hint = _did_you_mean(missing[1].strip('\'"'), layout_parts, 'layout parts')
                else:
                    hint = 'layout parts: ' + ', '.join(layout_parts)
                self.error(path, str(exc), hint)
            self._collect_encoding_lookups(path, ast.parse(expr.strip(), mode='eval'))

    def _collect_encoding_lookups(self, path, tree):
        parents = {child: node for node in ast.walk(tree) for child in ast.iter_child_nodes(node)}
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
                    and node.value.id == 'encoding'):
                continue
            key = node.slice
            if isinstance(key, ast.Name) and key.id == 'target':
                keys = None  # every target
            elif isinstance(key, ast.Constant) and isinstance(key.value, str):
                keys = [key.value]
            else:
                continue
            outer = parents.get(node)
            barcode = None
            if isinstance(outer, ast.Subscript) and outer.value is node and isinstance(outer.slice, ast.Constant) \
                    and isinstance(outer.slice.value, str):
                barcode = outer.slice.value
            self.lookups.append((path, keys, barcode))

    def _check_cycles(self):
        from .gen.probe import DAG
        try:
            dag = DAG()
            dag.from_config({'probes': self.protocol['probes']})
        except Exception as exc:
            self.error(('probes',), f'the probe graph cannot be built: {exc}',
                       'each probe needs template + parts; each part needs expr or template + parts')
            return
        _, cycles = _topological_order(dag.nodes)
        reported = set()
        for cycle in cycles:
            names = [node.name for node in cycle]
            if frozenset(names) in reported:
                continue
            reported.add(frozenset(names))
            at = next((name for name in names if self.nodes.get(name, ((),))[0][-1:] == ('expr',)), names[0])
            self.error(self.nodes.get(at, (('probes',),))[0], 'dependency cycle: ' + ' -> '.join(names + names[:1]),
                       'an expression cannot use itself or a probe that contains it; build it from target_region, '
                       'encoding or other probes')
        self.dag_ok = not cycles

    # -- encoding --------------------------------------------------------
    def check_encoding(self):
        if 'encoding' not in self.protocol:
            return
        encoding = self.protocol['encoding']
        if encoding is not None and not isinstance(encoding, dict):
            return  # reported by check_top_level
        encoding = encoding or {}
        needed = {}  # target -> ({barcode: None}, [expression paths])
        for path, keys, barcode in self.lookups:
            for key in (self.target_names if keys is None else keys):
                barcodes, users = needed.setdefault(key, ({}, []))
                if barcode is not None:
                    barcodes[barcode] = None
                users.append(format_path(path))
        for key, (barcodes, users) in needed.items():
            used_by = ', '.join(dict.fromkeys(users))
            if key not in encoding:
                entry = f'{key}:\n' + ''.join(f'  {b}: <sequence>\n' for b in barcodes or ['BC1'])
                close = _closest(key, encoding)
                prefix = f"did you mean '{close}'? Otherwise add" if close else 'add'
                self.error(('encoding',), f'no encoding entry for target {key!r} (used by {used_by})',
                           f'{prefix} under encoding:\n{entry.rstrip()}')
                continue
            entry = encoding[key]
            if barcodes and not isinstance(entry, dict):
                self.error(('encoding', key), f'encoding of {key!r} must map barcode names to sequences',
                           'e.g.\n' + '\n'.join(f'{b}: <sequence>' for b in barcodes))
                continue
            for barcode in barcodes:
                if barcode not in entry:
                    hint = _did_you_mean(barcode, entry, 'barcodes')
                    self.error(('encoding', key), f'target {key!r} has no barcode {barcode!r} (used by {used_by})',
                               (hint + '; or ' if hint else '') + f'add {barcode}: <sequence>')
                elif isinstance(entry[barcode], list):
                    # A list of sequences, indexed by the expression (e.g. with probe_index).
                    for i, value in enumerate(entry[barcode]):
                        self._check_barcode_sequence(('encoding', key, barcode, i), f'{barcode}[{i}]', key, value)
                    if not entry[barcode]:
                        self.error(('encoding', key, barcode), f'barcode list {barcode!r} of {key!r} is empty',
                                   f'list the sequences, e.g. {barcode}: [<sequence>, <sequence>]')
                else:
                    self._check_barcode_sequence(('encoding', key, barcode), barcode, key, entry[barcode])
        if any(keys is None for _, keys, _ in self.lookups):
            for key in encoding:
                if key not in self.target_names:
                    self.warning(('encoding', key), f'encoding entry {key!r} is not a target and is never used',
                                 _did_you_mean(key, self.target_names, 'targets'))

    def _check_barcode_sequence(self, path, barcode, key, value):
        if not isinstance(value, str) or not value.strip():
            self.error(path, f'barcode {barcode!r} of {key!r} is empty', f'set {barcode} to a nucleotide sequence')
        elif set(value.upper()) - IUPAC:
            self.warning(path, f'barcode {barcode!r} of {key!r} is not a nucleotide sequence: {value!r}',
                         'use only A, C, G, T/U and IUPAC codes such as N')

    # -- hybridization conditions ----------------------------------------
    def check_hybridization(self):
        """Check the protocol block and each attribute override on its own, so every one is reported."""
        blocks = []
        if 'hybridization' in self.protocol:
            blocks.append((('hybridization',), {'hybridization': self.protocol['hybridization']}))
        attributes = self.protocol.get('attributes')
        for name, attr in (attributes.items() if isinstance(attributes, dict) else []):
            if isinstance(attr, dict) and 'tm_method' in attr:
                try:
                    method = tm_method(attr, f'attributes.{name}')
                    try:
                        ignored = unused_fields(method, resolve_conditions(self.protocol, name))
                    except ValueError:
                        ignored = []  # the hybridization block itself is reported below
                    if ignored:
                        self.warning(('attributes', name, 'tm_method'),
                                     f'tm_method {method} does not use {", ".join(ignored)}, declared for this attribute',
                                     'remove the field(s) or choose another tm_method')
                except ValueError as exc:
                    hint = (f'tm_method: {TM_METHODS[1]}   # or {", ".join(TM_METHODS[2:])}, {TM_METHODS[0]} (the default)'
                            if 'must be one of' in str(exc) else 'remove it')
                    self.error(('attributes', name, 'tm_method'), str(exc), hint)
                attr = {key: value for key, value in attr.items() if key != 'tm_method'}
            if isinstance(attr, dict) and 'hybridization' in attr:
                blocks.append((('attributes', name, 'hybridization'), {'attributes': {name: attr}}))
        for path, protocol in blocks:
            try:
                validate_protocol_hybridization(protocol)
            except ValueError as exc:
                hint = ('remove it; conditions apply to Tm and fold_score attributes or to the whole protocol '
                        '(top-level hybridization)') if 'only supported' in str(exc) else HYBRIDIZATION_HINT
                self.error(path, str(exc), hint)

    def check_default_tm_filters(self, protocol):
        """Generated Tm filters start at the hybridization temperature, or at 37 °C when none is declared."""
        filters = (protocol.get('post_process') or {}).get('filters')
        assumed = {}
        for name in (filters if isinstance(filters, dict) else {}):
            if not name.endswith('_tm'):
                continue
            try:
                low, high, source = default_tm_window(protocol, name)
            except ValueError:
                return  # invalid conditions are reported by check_hybridization
            if source == 'default':
                assumed[name] = f'{low:g}-{high:g} °C'
        if assumed:
            self.warning(('hybridization',), 'salt, probe or formamide conditions are set but temperature_c is not, '
                         f'so the generated Tm filters ({", ".join(assumed)}) assume the default 37 °C and keep '
                         f'{next(iter(assumed.values()))}',
                         'add temperature_c: <your hybridization temperature>, or write post_process.filters')

    # -- runtime normalization -------------------------------------------
    def normalize(self):
        """Run the normalization of `uprobe run`; return (protocol, auto-generated sections)."""
        p = self.protocol
        if any(key not in p for key in REQUIRED_KEYS):
            return p, set()
        from .cli import _is_blank_step, _validate_and_normalize_protocol
        generated = {key for key in ('attributes', 'summary') if not p.get(key)}
        rna = self.region().get('source', 'genome') != 'genome'
        post_process = p.get('post_process')
        if isinstance(post_process, dict):
            post_process = {k: v for k, v in post_process.items()
                            if not _is_blank_step(v) and not (rna and k in DNA_ONLY_STEPS)}
        if not post_process:
            generated.add('post_process')
        logger = logging.getLogger('uprobe.core.cli')
        level = logger.level
        logger.setLevel(logging.WARNING)
        try:
            # As `uprobe run`, which picks the aligner of a generated off-target attribute by genome.
            genome = self.genomes.get(p['genome']) if isinstance(self.genomes, dict) \
                and isinstance(p.get('genome'), str) else None
            return _validate_and_normalize_protocol(copy.deepcopy(p), genome), generated
        except Exception as exc:
            if not any(f.severity == 'error' for f in self.findings):  # otherwise already explained
                self.error((), f'`uprobe run` rejects this protocol: {exc}', 'fix the section named in the message')
            return p, set()
        finally:
            logger.setLevel(level)

    # -- attributes ------------------------------------------------------
    def check_attributes(self, protocol, auto):
        attributes = protocol.get('attributes')
        self.attributes = attributes if isinstance(attributes, dict) else {}
        if auto or not self.attributes:
            return  # generated attributes are only checked against the genome config
        columns = self.sequence_columns()
        index_nodes = probe_index_nodes(protocol.get('probes'))
        for name, attr in self.attributes.items():
            path = ('attributes', name)
            if not isinstance(attr, dict):
                self.error(path, f'attribute {name!r} must be a mapping with `target` and `type`',
                           f'e.g.\n{name}:\n  target: target_region\n  type: gc_content')
                continue
            for option in attr:
                if option not in ATTRIBUTE_KEYS:
                    self.warning(path + (option,), f'attribute option {option!r} is ignored',
                                 _did_you_mean(option, ATTRIBUTE_KEYS, 'known options'))
            kind = attr.get('type')
            kind = str(kind).lower() if kind else None
            if kind is None:
                self.error(path + ('type',), 'missing attribute `type`', 'valid types: ' + ', '.join(ATTRIBUTE_TYPES))
            elif kind not in ATTRIBUTE_TYPES:
                hint = (f"did you mean '{ATTRIBUTE_ALIASES[kind]}'? ({kind} was renamed)" if kind in ATTRIBUTE_ALIASES
                        else _did_you_mean(kind, ATTRIBUTE_TYPES, 'valid types'))
                self.error(path + ('type',), f'unknown attribute type {attr["type"]!r}', hint)
            else:
                aligners = ATTRIBUTE_TYPES[kind]
                if aligners and attr.get('aligner') not in aligners:
                    got = f', got {attr["aligner"]!r}' if attr.get('aligner') else ''
                    self.error(path + ('aligner',), f'{kind} requires aligner {" or ".join(aligners)}{got}',
                               f'aligner: {aligners[0]}')
                elif kind == 'mapped_genes':
                    self._check_mapped_genes(path, attr)
                elif kind == 'kmer_count' and attr.get('strand', 'both') not in STRANDS:
                    self.error(path + ('strand',), f'strand must be one of {", ".join(STRANDS)}, got '
                               f'{attr["strand"]!r}', KMER_STRAND_HINT)
                if kind == 'blocks' and self.source == 'genome' and not self.direct_targets:
                    self.warning(path + ('type',), 'blocks needs a start column, which source: genome does not '
                                 'produce; the attribute is skipped', 'remove it, or use source mRNA, exon, CDS or UTR')
            if kind == 'difference':
                self._check_difference(path, attr, list(self.attributes)[:list(self.attributes).index(name)])
                continue
            target = attr.get('target')
            if target is None:
                self.error(path + ('target',), 'missing attribute `target`',
                           'target: target_region, or a probe/part such as probe_1.part1')
            elif not isinstance(target, str) or (target not in columns and target.replace(':', '.') not in columns):
                self.error(path + ('target',), f'target {target!r} is not target_region, a probe/part or a layout '
                           'part; the attribute would be skipped',
                           _did_you_mean(str(target).replace(':', '.'), columns, 'valid targets'))
            elif target.replace(':', '.') in index_nodes:
                self.warning(path + ('target',), f'{target!r} depends on probe_index, which numbers the final '
                             'probes after post-processing: filters and sorts on this attribute see each candidate '
                             'numbered in extraction order; the final table recomputes it',
                             'compute it on the parts that do not read probe_index, or accept the approximation')
            if name in columns or name in TARGET_COLUMNS:
                self.warning(path, f'attribute name {name!r} overwrites the column of the same name',
                             f'rename it, e.g. {name}_{kind or "value"}')

    def _check_difference(self, path, attr, earlier):
        of = attr.get('of')
        if not (isinstance(of, list) and len(of) == 2 and all(isinstance(name, str) for name in of)):
            self.error(path + ('of',), f'difference needs `of`: two attribute names, got {of!r}',
                       'of: [padlock_tm, primer_tm]   # |padlock_tm - primer_tm|')
            return
        for name in of:
            if name not in earlier:
                self.error(path + ('of',), f'{name!r} is not an attribute defined before this one',
                           _did_you_mean(name, earlier, 'attributes defined before it'))

    def _check_mapped_genes(self, path, attr):
        hints = {'strand': 'strand: sense   # both, sense (for target_region) or antisense (for a probe binding part)',
                 'max_mismatches': f'max_mismatches: 2   # substitutions, 0-{MAX_MISMATCHES_LIMIT}',
                 'exempt_copies': COPIES_SNIPPET}
        for option, hint in hints.items():
            if option in attr:
                try:
                    mapped_genes_settings({'aligner': attr['aligner'], option: attr[option]})
                except ValueError as exc:
                    self.error(path + (option,), str(exc), hint)
        if attr['aligner'] == 'bowtie2' and 'max_mismatches' in attr:
            self.warning(path + ('max_mismatches',), 'max_mismatches only applies to aligner: bowtie and is ignored',
                         'aligner: bowtie   # every site within max_mismatches substitutions (Bowtie 1)')
        step = self.region().get('step')
        if attr.get('exempt_copies') and self.lengths and type(step) is int and step > self.lengths[1]:
            self.warning(path + ('exempt_copies',), f'extracts.target_region.step {step} leaves gaps between the '
                         f'windows of at most {self.lengths[1]} nt, and exempt_copies reads the target\'s sequence '
                         'from the windows: k-mers across the gaps are unknown, so copies are judged on part of the '
                         'target or not found', f'step: {self.lengths[1]}   # or less, so the windows cover the target')

    # -- gene_expression ---------------------------------------------------
    def check_gene_expression(self, protocol):
        from .gene_expression import identifier_level, read_table, resolve_columns, validate_gene_expression
        if self.protocol.get('gene_expression') is None:
            return
        path = ('gene_expression',)
        try:
            settings = validate_gene_expression(self.protocol['gene_expression'])
        except ValueError as exc:
            self.error(path, str(exc), 'e.g.\n' + GENE_EXPRESSION_SNIPPET)
            return
        attributes = protocol.get('attributes')
        if not any(isinstance(attr, dict) and str(attr.get('type', '')).lower() == 'mapped_genes'
                   for attr in (attributes.values() if isinstance(attributes, dict) else [])):
            self.warning(path, 'gene_expression only changes how mapped_genes attributes count off-target genes, '
                         'and there is none', 'add e.g.\ntarget_mapped_genes:\n  target: target_region\n'
                         '  type: mapped_genes\n  aligner: bowtie')
        table = Path(settings['table'])
        if not table.is_file():
            self.error(path + ('table',), f'expression table not found: {table}',
                       'fix the path; relative paths are resolved from the directory where uprobe runs')
            return
        try:
            df = read_table(table)
        except Exception as exc:
            self.error(path + ('table',), f'cannot read {table}: {exc}',
                       'use a tab-separated table (comma-separated when the name contains .csv) with a header row')
            return
        try:
            id_columns, _ = resolve_columns(df, settings)
        except ValueError as exc:
            key = next((key for key in ('value_column', 'transcript_column') if key in str(exc)), 'gene_column')
            self.error(path + (key,), f'{table}: {exc}', f'{key}: one of {", ".join(map(str, df.columns))}')
            return
        # Transcript ids match no gene, and gene ids no transcript.
        level = settings['level']
        other = 'transcript' if level == 'gene' else 'gene'
        if all(identifier_level(df[column].tolist()) == other for column in id_columns):
            example = next((str(v) for v in df[id_columns[0]] if isinstance(v, str) and v.strip()), '')
            columns = ', '.join(id_columns)
            replaced = self.protocol['gene_expression'].get(f'{level}_column') is not None
            self.error(path + ('level',), f'{table}: level is {level}, but column {columns} holds {other} ids '
                       f'(e.g. {example}), which match no {level} of the annotation',
                       f'level: {other}' + (f'\n{other}_column: {id_columns[0]}   # in place of {level}_column'
                                            if replaced else ''))

    # -- post_process ----------------------------------------------------
    def check_post_process(self, protocol):
        from .cli import _is_blank_step
        raw = self.protocol.get('post_process')
        rna = self.source not in (None, 'genome')
        if isinstance(raw, dict) and rna:
            for key in DNA_ONLY_STEPS:
                if raw.get(key):
                    self.warning(('post_process', key), f'{key} only applies to source: genome and is dropped '
                                 f'for source: {self.source}', f'remove post_process.{key}')
        for key, value in (raw.items() if isinstance(raw, dict) else []):
            if isinstance(value, dict) and value and _is_blank_step(value):
                self.warning(('post_process', key), f'every option of {key} is blank, so `uprobe run` skips it',
                             'fill in its options (e.g. remove_overlap: {location_interval: 0}), or delete it')
        post_process = protocol.get('post_process')
        if not isinstance(post_process, dict) or not post_process:
            return
        # As in _validate_and_normalize_protocol, when normalization was skipped.
        post_process = {k: v for k, v in post_process.items()
                        if not _is_blank_step(v) and not (rna and k in DNA_ONLY_STEPS)}
        columns = self.known_columns(protocol)
        for key in post_process:
            if key not in POST_PROCESS_KEYS:
                # A misspelled step silently disables that step.
                self.error(('post_process', key), f'unknown post_process step {key!r} is ignored',
                           _did_you_mean(key, POST_PROCESS_KEYS, 'known steps'))
        self._check_filters(post_process.get('filters'), columns)
        self._check_sorts(post_process.get('sorts'), columns)
        self._check_remove_overlap(post_process.get('remove_overlap'))
        self._check_equal_space(post_process.get('equal_space'), post_process.get('remove_overlap'))
        self._check_avoid_otp(post_process.get('avoid_otp'))
        self._check_off_target_pool(post_process.get('off_target_pool'), post_process.get('filters'))

    def _check_filters(self, filters, columns):
        path = ('post_process', 'filters')
        if filters is None:
            return
        if not isinstance(filters, dict):
            self.error(path, 'filters must be a mapping of filter name -> {condition: ...}',
                       'e.g.\n' + SNIPPETS['post_process'])
            return
        for name, config in filters.items():
            fpath = path + (name,)
            if not isinstance(config, dict):
                self.error(fpath, f'filter {name!r} must be a mapping with `condition`',
                           f'e.g.\n{name}:\n  condition: {name} >= 50 & {name} <= 65')
            elif config.get('type') == 'sequence_pattern':
                self._check_sequence_filter(fpath, config)
            elif 'type' in config:
                self.error(fpath + ('type',), f'unknown filter type {config["type"]!r}',
                           'use type: sequence_pattern for sequence exclusions, or omit type for a numeric condition')
            elif 'condition' not in config:
                self.error(fpath, f'filter {name!r} has no `condition` and is skipped',
                           f'condition: {name} >= 50 & {name} <= 65')
            elif isinstance(config['condition'], list):
                for i, condition in enumerate(config['condition']):
                    self._check_condition(fpath + ('condition', i), condition, columns)
            else:
                self._check_condition(fpath + ('condition',), config['condition'], columns)

    def _check_condition(self, path, condition, columns):
        if not isinstance(condition, str):
            self.error(path, f'condition must be a string, got {condition!r}', CONDITION_HINT)
            return
        try:
            parsed = parse_condition(condition)
        except ValueError as exc:
            self.error(path, str(exc), CONDITION_HINT)
            return
        for (column, operator, value), text in zip(parsed, condition.split('&')):
            text = text.strip()
            if operator not in OPERATORS:
                fix = OPERATOR_FIXES.get(operator)
                self.error(path, f'unsupported operator {operator!r} in {text!r}',
                           f"use '{fix}'" if fix else 'use one of ' + ', '.join(OPERATORS))
                continue
            rest = re.sub(r'\s+', '', text)[len(column) + len(operator):]
            try:
                number = float(rest)
            except ValueError:
                self.error(path, f'cannot read {text!r}: only "<column> <op> <number>" is understood and the rest '
                           'is ignored', CONDITION_HINT + '; there is no "or": use separate filters')
            else:
                if number != value:
                    self.error(path, f'{text!r} is read as {column} {operator} {value}',
                               'write the number in plain decimal notation, e.g. 0.00001')
            if column not in columns:
                self.error(path, f'condition uses unknown column {column!r}; this comparison is silently skipped',
                           _did_you_mean(column, columns, 'attributes'))

    def _check_sequence_filter(self, path, config):
        import pandas as pd
        from .process.sequence_filter import apply_sequence_filter
        target = config.get('target')
        columns = self.sequence_columns()
        if isinstance(target, str) and target and target not in columns and target.replace(':', '.') not in columns:
            self.error(path + ('target',), f'sequence filter target {target!r} is not a sequence column',
                       _did_you_mean(target.replace(':', '.'), columns, 'sequence columns'))
        column = target if isinstance(target, str) and target else 'target_region'
        try:
            apply_sequence_filter(pd.DataFrame({column: ['ACGT']}), config)
        except ValueError as exc:
            self.error(path, str(exc), 'e.g.\ntype: sequence_pattern\ntarget: probe_1\n'
                       'exclude_patterns: ["^[GC]{2}", "AAAA"]')

    def _check_sorts(self, sorts, columns):
        path = ('post_process', 'sorts')
        if sorts is None:
            return
        if not isinstance(sorts, dict):
            self.error(path, 'sorts must be a mapping with is_ascending / is_descending lists',
                       'sorts:\n  is_ascending: [probe_1_fold]\n  is_descending: [probe_1_tm]')
            return
        for key in sorts:
            if key not in SORT_KEYS:
                self.error(path + (key,), f'unknown sort key {key!r} is ignored', _did_you_mean(key, SORT_KEYS))
        seen = {}
        for key in SORT_KEYS:
            fields = sorts.get(key)
            if fields is None:
                continue
            if not isinstance(fields, list):
                self.error(path + (key,), f'{key} must be a list of column names', f'{key}: [{fields}]')
                continue
            for i, name in enumerate(fields):
                if not isinstance(name, str) or name not in columns:
                    self.error(path + (key, i), f'unknown sort column {name!r}; it is dropped and shifts the '
                               'direction of the columns after it', _did_you_mean(name, columns, 'attributes'))
                elif seen.setdefault(name, key) != key:
                    self.warning(path + (key, i), f'{name!r} is in both is_ascending and is_descending',
                                 'keep it in one list')
                elif key == 'is_descending' and self._attribute_type(name) in LOWER_IS_BETTER:
                    worst = LOWER_IS_BETTER[self._attribute_type(name)]
                    self.warning(path + (key, i), f'{name!r} is sorted descending, which ranks the {worst} '
                                 'candidates first; lower is better', f'move {name} to is_ascending')

    def _attribute_type(self, name):
        attr = self.attributes.get(name)
        kind = attr.get('type') if isinstance(attr, dict) else None
        return ATTRIBUTE_ALIASES.get(kind.lower(), kind.lower()) if isinstance(kind, str) else None

    def _check_remove_overlap(self, config):
        path = ('post_process', 'remove_overlap')
        if config is None:
            return
        if not isinstance(config, dict):
            self.error(path, f'remove_overlap must be a mapping, got {config!r}',
                       'remove_overlap:\n  location_interval: 0')
            return
        for key in config:
            if key not in ('location_interval', 'prefer'):
                self.warning(path + (key,), f'unknown option {key!r} is ignored',
                             'the options are location_interval and prefer')
        prefer = config.get('prefer')
        if prefer is not None and (not isinstance(prefer, str) or prefer not in self.attributes):
            self.error(path + ('prefer',), f'prefer must name an attribute (lower is better), got {prefer!r}',
                       _did_you_mean(str(prefer), list(self.attributes), 'attributes'))
        interval = config.get('location_interval', 0)
        if type(interval) is not int:
            self.error(path + ('location_interval',), f'location_interval must be an integer, got {interval!r}',
                       'location_interval: 0   # minimum gap in bases between kept probes')

    def _check_positive_int(self, path, value):
        if type(value) is not int or value <= 0:
            self.error(path, f'number_desired must be a positive integer, got {value!r}', 'number_desired: 50')

    def _check_equal_space(self, config, remove_overlap):
        path = ('post_process', 'equal_space')
        if config is None:
            return
        if not isinstance(config, dict):
            self.error(path, 'equal_space must be a mapping',
                       'equal_space:\n  number_desired: 50   # or per target: {TARGET: {number_desired: 50}}')
            return
        for key, value in config.items():
            if key == 'number_desired':
                self._check_positive_int(path + (key,), value)
                continue
            if self.target_names and key not in self.target_names:
                self.warning(path + (key,), f'{key!r} is not a target; this entry is ignored',
                             _did_you_mean(key, self.target_names, 'targets'))
            if not isinstance(value, dict) or 'number_desired' not in value:
                self.error(path + (key,), f'equal_space entry {key!r} needs number_desired',
                           f'{key}:\n  number_desired: 50')
            else:
                self._check_positive_int(path + (key, 'number_desired'), value['number_desired'])
                for option in value:
                    if option != 'number_desired':
                        self.warning(path + (key, option), f'unknown option {option!r} is ignored',
                                     'the only option is number_desired')
        if 'number_desired' not in config:
            unset = [name for name in self.target_names if name not in config]
            if unset:
                self.warning(path, f'target(s) {", ".join(map(repr, unset))} have no number_desired and keep '
                             'all their probes', 'add a global number_desired, or an entry per target')
        region = self.region()
        try:
            _, _, step = normalize_sampling(region.get('length'), region.get('step'), region.get('overlap'))
        except ValueError:
            step = None  # reported under extracts
        if remove_overlap is None and step is not None and self.lengths and step < self.lengths[0]:
            self.warning(path, 'without remove_overlap, overlapping candidate windows can be selected',
                         'add remove_overlap:\n  location_interval: 2   # minimum gap in bases between probes')

    def _check_off_target_pool(self, config, filters):
        from .process.off_target_pool import OPTIONS, pool_settings
        path = ('post_process', 'off_target_pool')
        if config is None:
            return
        snippet = ('off_target_pool:\n  attribute: target_mapped_genes   # a mapped_genes attribute\n'
                   '  max_probes_per_gene: 2   # probes of one target that may hit the same off-target gene')
        if not isinstance(config, dict):
            self.error(path, 'off_target_pool must be a mapping', snippet)
            return
        for key in config:
            if key not in OPTIONS:
                self.error(path + (key,), f'unknown option {key!r} of off_target_pool',
                           _did_you_mean(key, OPTIONS, 'options'))
        limit = config.get('max_probes_per_gene')
        if type(limit) is not int or limit < 0:
            self.error(path + ('max_probes_per_gene',), f'max_probes_per_gene must be an integer >= 0, got {limit!r}',
                       'max_probes_per_gene: 2   # 0 keeps only probes without off-target genes')
        try:
            names = pool_settings({'attribute': config.get('attribute'), 'max_probes_per_gene': 0})[0]
        except ValueError as exc:
            self.error(path + ('attribute',), str(exc), snippet)
            return
        mapped = [name for name, attr in self.attributes.items()
                  if isinstance(attr, dict) and str(attr.get('type', '')).lower() == 'mapped_genes']
        for name in names:
            if name not in mapped:
                self.error(path + ('attribute',), f'{name!r} is not a mapped_genes attribute, whose off-target '
                           'genes the limit counts', _did_you_mean(name, mapped, 'mapped_genes attributes')
                           or 'add one, e.g.\ntarget_mapped_genes:\n  target: target_region\n  type: mapped_genes\n'
                              '  aligner: bowtie')
            elif self._filters_out_off_targets(name, filters):
                self.warning(path, f'the filter {name} <= 1 already removes every probe with an off-target gene '
                             '(filters run first), so off_target_pool removes nothing',
                             f'relax the filter on {name} (e.g. {name} <= 4) or remove it')

    @staticmethod
    def _filters_out_off_targets(name, filters):
        """True if a filter keeps only rows where mapped_genes attribute `name` is at most 1."""
        for config in (filters.values() if isinstance(filters, dict) else []):
            conditions = config.get('condition') if isinstance(config, dict) else None
            for condition in conditions if isinstance(conditions, list) else [conditions]:
                try:
                    parsed = parse_condition(condition) if isinstance(condition, str) else []
                except ValueError:
                    continue
                if any(column == name and ((op in ('<=', '==') and value <= 1) or (op == '<' and value <= 2))
                       for column, op, value in parsed):
                    return True
        return False

    def _check_avoid_otp(self, config):
        from .process.otp import parse_region
        path = ('post_process', 'avoid_otp')
        if config is None:
            return
        if not isinstance(config, dict):
            self.error(path, 'avoid_otp must be a mapping of target -> settings', OTP_SNIPPET)
            return
        if not any(isinstance(a, dict) and str(a.get('type', '')).lower() == 'mapped_sites'
                   for a in self.attributes.values()):
            self.error(path, 'avoid_otp reads the sites of a mapped_sites attribute and is skipped without one',
                       'add e.g.\nprobe_1_sites:\n  target: probe_1\n  type: mapped_sites\n  aligner: bowtie2')
        for key, value in config.items():
            tpath = path + (key,)
            if self.target_names and key not in self.target_names:
                self.warning(tpath, f'{key!r} is not a target; this entry is ignored',
                             _did_you_mean(key, self.target_names, 'targets'))
            if not isinstance(value, dict):
                self.error(tpath, f'avoid_otp entry {key!r} must be a mapping', OTP_SNIPPET)
                continue
            regions = value.get('target_regions', [])
            if not isinstance(regions, list):
                self.error(tpath + ('target_regions',), 'target_regions must be a list of chrom:start-end regions',
                           'target_regions: [chr1:1000-2000]')
                regions = []
            for i, region in enumerate(regions):
                try:
                    parse_region(str(region))
                except ValueError:
                    self.error(tpath + ('target_regions', i), f'invalid region {region!r}',
                               'write chrom:start-end, e.g. chr1:100000-101000')
            try:
                float(value.get('density_thresh', 1e-4))
            except (TypeError, ValueError):
                self.error(tpath + ('density_thresh',), f'density_thresh must be a number, got '
                           f'{value["density_thresh"]!r}', 'density_thresh: 0.0001')
            try:
                search_range = tuple(map(float, value.get('search_range', [-1000000, 1000000])))
                if len(search_range) != 2:
                    raise ValueError
            except (TypeError, ValueError):
                self.error(tpath + ('search_range',), 'search_range must be two numbers',
                           'search_range: [-1000000, 1000000]')

    # -- summary ---------------------------------------------------------
    def check_summary(self, protocol):
        summary = protocol.get('summary')
        if not isinstance(summary, dict) or not summary:
            return
        for key in summary:
            if key not in SUMMARY_KEYS:
                self.warning(('summary', key), f'unknown summary key {key!r} is ignored',
                             _did_you_mean(key, SUMMARY_KEYS, 'known keys'))
        if not summary.get('report_name'):
            self.warning(('summary', 'report_name'), 'report_name is empty, so no HTML report is generated',
                         'report_name: ' + ('dna_report' if self.source == 'genome' else 'rna_report'))
        attributes = summary.get('attributes')
        if attributes is None:
            return
        if isinstance(attributes, dict):
            items = [(('summary', 'attributes', name), name, config) for name, config in attributes.items()]
        elif isinstance(attributes, list):
            items = [(('summary', 'attributes', i), name, None) for i, name in enumerate(attributes)]
        else:
            self.error(('summary', 'attributes'), 'summary.attributes must be a list of attribute names',
                       'attributes:\n  - probe_1_tm')
            return
        columns = self.known_columns(protocol)
        for path, name, config in items:
            if not isinstance(name, str) or name not in columns:
                self.error(path, f'summary attribute {name!r} is not an attribute or column, so the report leaves '
                           'it out', _did_you_mean(name, columns, 'attributes'))
            stats = config.get('stats') if isinstance(config, dict) else None
            for stat in stats if isinstance(stats, list) else []:
                if stat not in SUMMARY_STATS:
                    self.warning(path + ('stats',), f'unknown statistic {stat!r} is ignored',
                                 _did_you_mean(stat, SUMMARY_STATS, 'valid statistics'))

    # -- genome config ---------------------------------------------------
    def check_genome(self, protocol, auto):
        genomes, doc = self.genomes, 'genomes'
        if not isinstance(genomes, dict) or not genomes:
            self.error((), 'the genome config must map genome names to {fasta, gtf, align_index}',
                       'e.g.\n' + GENOME_SNIPPET, doc)
            return
        name = self.protocol.get('genome')
        if not isinstance(name, str) or not name:
            return  # reported as a protocol error
        if name not in genomes:
            self.error(('genome',), f'genome {name!r} is not defined in the genome config',
                       _did_you_mean(name, genomes, 'available genomes'))
            return
        entry = genomes[name]
        if not isinstance(entry, dict):
            self.error((name,), f'genome {name!r} must be a mapping', 'e.g.\n' + GENOME_SNIPPET, doc)
            return
        files_ok = True
        for key, ext in (('fasta', 'fa'), ('gtf', 'gtf')):
            value = entry.get(key)
            if not value:
                self.error((name, key) if key in entry else (name,), f'genome {name!r} has no {key} path',
                           f'{key}: /path/to/genome.{ext}', doc)
                files_ok = False
            elif not Path(str(value)).exists():
                self.error((name, key), f'{key} file not found: {value}',
                           'fix the path; relative paths are resolved from the directory where uprobe runs', doc)
                files_ok = False
        align = entry.get('align_index')
        if align is None:
            self.error((name,), f'genome {name!r} has no align_index; `uprobe run` reads it to build indexes',
                       'align_index:\n  - bowtie2', doc)
            align = []
        elif not isinstance(align, list):
            self.error((name, 'align_index'), 'align_index must be a list', 'align_index:\n  - bowtie2', doc)
            align = []
        for i, aligner in enumerate(align):
            if aligner not in ALIGN_INDEXES:
                self.error((name, 'align_index', i), f'unknown aligner {aligner!r}',
                           _did_you_mean(aligner, ALIGN_INDEXES, 'supported aligners'), doc)
            else:
                missing = [exe for exe in ALIGN_INDEXES[aligner] if not shutil.which(exe)]
                if missing:
                    self.warning((name, 'align_index', i), f'{" and ".join(missing)} not found on PATH; indexing and '
                                 f'aligning with {aligner} will fail',
                                 f'install it, e.g. conda install -c bioconda {ALIGN_PACKAGES[aligner]}', doc)
        attributes = protocol.get('attributes')
        needs = {'bowtie2': [], 'bowtie': [], 'jellyfish': []}
        for attr_name, attr in (attributes.items() if isinstance(attributes, dict) else []):
            aligners = ATTRIBUTE_TYPES.get(str(attr.get('type', '')).lower()) if isinstance(attr, dict) else None
            if aligners:
                needs[attr.get('aligner') if attr.get('aligner') in aligners else aligners[0]].append(attr_name)
        enabled = {'bowtie2': 'bowtie2' in align, 'bowtie': 'bowtie' in align,
                   'jellyfish': entry.get('jellyfish', False) is True}
        fixes = {'bowtie2': f'add bowtie2 to align_index of {name!r} in the genome config',
                 'bowtie': f'add bowtie to align_index of {name!r} in the genome config '
                           f'(`uprobe build-index` then builds its transcriptome index; Bowtie 1: {BOWTIE_INSTALL})',
                 'jellyfish': f'set jellyfish: true for {name!r} in the genome config'}
        for required, names in needs.items():
            if not names:
                continue
            if not enabled[required]:
                if auto:
                    self.error(('attributes',), f'the attributes generated because `attributes` is empty '
                               f'({", ".join(names)}) need {required}, which genome {name!r} does not enable',
                               fixes[required] + ', or define attributes explicitly')
                else:
                    for attr_name in names:
                        self.error(('attributes', attr_name), f'attribute {attr_name!r} needs {required}, which '
                                   f'genome {name!r} does not enable', fixes[required])
            elif required == 'jellyfish' and not shutil.which('jellyfish'):
                self.warning((name, 'jellyfish'), 'jellyfish is not on PATH; kmer_count attributes will fail',
                             'install it, e.g. conda install -c bioconda jellyfish', doc)
        if self.lookup_targets and files_ok:
            self.check_targets_in_genome(entry)

    def check_targets_in_genome(self, entry):
        if self.source == 'genome':
            try:
                lengths = _fasta_lengths(entry['fasta'])
            except Exception as exc:
                self.error(('genome',), f'cannot read the FASTA index: {exc}', 'check the fasta file')
                return
            for target in self.genome_targets:
                match = REGION_RE.match(str(target))
                if not match:
                    continue
                path = self.target_paths.get(target, ('targets',))
                chrom, end = match[1], int(match[3])
                if chrom not in lengths:
                    self.error(path, f'chromosome {chrom!r} of target {target!r} is not in the FASTA',
                               _did_you_mean(chrom, lengths, 'FASTA sequences'))
                elif end > lengths[chrom]:
                    self.error(path, f'target {target!r} ends after the end of {chrom} ({lengths[chrom]} bp)',
                               f'use an end <= {lengths[chrom]}')
            if self.region().get('masked') == 'drop':
                self.check_masked_genome(entry['fasta'], lengths)
        elif self.source is not None and self.genome_targets:
            from .gen.fun import read_gtf
            try:
                # The same gene_name lookup as gen.fun.validate_targets.
                transcripts = self.region().get('transcripts', 'all')
                gtf = read_gtf(entry['gtf'], extract_fields=['gene_name'] if transcripts == 'all'
                               else ['gene_name', 'transcript_id'])
                genes = set(gtf['gene_name'].dropna().unique())
            except Exception as exc:
                self.error(('genome',), f'cannot read the GTF: {exc}', 'check the gtf file')
                return
            coding = set(gtf.loc[gtf['type'] == 'CDS', 'gene_name'].dropna())
            self.coding_targets = [target for target in self.genome_targets if target in coding]
            self.check_transcripts_in_genome(gtf, [target for target in self.genome_targets if target in genes])
            for target in self.genome_targets:
                if target not in genes:
                    close = _closest(target, genes)
                    self.error(self.target_paths.get(target, ('targets',)),
                               f'target {target!r} is not a gene_name in the GTF',
                               f"did you mean '{close}'?" if close else
                               'check the gene symbol (GTF gene_name/gene attributes are matched, not gene IDs)')

    def check_masked_genome(self, fasta, lengths):
        """Warn when masked: drop reads targets without a soft-masked base: the FASTA is probably not repeat-masked."""
        from pyfaidx import Fasta
        regions = [m for m in (REGION_RE.match(str(t)) for t in self.genome_targets)
                   if m and m[1] in lengths and int(m[3]) <= lengths[m[1]]]
        if not regions:
            return
        with Fasta(str(fasta)) as fa:
            if any(not fa[m[1]][int(m[2]) - 1:int(m[3])].seq.isupper() for m in regions):
                return
        self.warning(('extracts', 'target_region', 'masked'), 'the targets have no soft-masked (lowercase) base: '
                     'the FASTA looks unmasked, so masked: drop skips only windows with N, not repeats',
                     'use a soft-masked FASTA, e.g. `uprobe genome fetch --species homo_sapiens --soft-masked` '
                     'or UCSC hg38.fa')

    def check_transcripts_in_genome(self, gtf, targets):
        """The transcripts option against the GTF: listed ids, and the tags canonical and mane read."""
        from .gen.fun import TranscriptSelectionError, select_transcripts, transcripts_setting
        try:
            setting = transcripts_setting(self.region().get('transcripts'))
        except ValueError:
            return  # reported by check_transcripts
        if setting == 'all' or not targets:
            return
        logger = logging.getLogger('uprobe.core.gen.fun')
        level = logger.level
        logger.setLevel(logging.WARNING)
        try:
            select_transcripts(gtf, targets, setting)
        except TranscriptSelectionError as exc:
            for problem in exc.problems:
                self.error(('extracts', 'target_region', 'transcripts'), problem, TRANSCRIPTS_HINT)
        finally:
            logger.setLevel(level)

    def check_exon_source(self):
        """Warn that source: exon is the coding sequence of protein-coding genes (see gen.fun.RNA_SOURCES)."""
        if self.source != 'exon' or not self.genome_targets:
            return
        path = ('extracts', 'target_region', 'source')
        fix = ('source: mRNA   # the whole mature mRNA: 5\' UTR, CDS and 3\' UTR\n'
               'or source: CDS to keep the coding sequence and say so')
        if self.coding_targets is None:  # the GTF was not read
            self.warning(path, 'source: exon extracts only the coding sequence (CDS features) of protein-coding '
                         'genes, without their UTRs; only genes without CDS keep all their exons', fix)
        elif self.coding_targets:
            shown = ', '.join(self.coding_targets[:10]) + (', ...' if len(self.coding_targets) > 10 else '')
            self.warning(path, f'{len(self.coding_targets)} target(s) are protein-coding ({shown}): source: exon '
                         'extracts only their coding sequence (CDS features), without the UTRs', fix)

    # -- probe construction on synthetic targets -------------------------
    def dry_run(self):
        # ExprProbe.build uses the safe expression evaluator, so building the
        # probes never runs arbitrary code.
        if self.has_errors('probes', 'encoding', 'extracts') or not self.dag_ok \
                or not self.target_names or self.lengths is None:
            return
        from .gen.probe import DAG
        targets = self.target_names[:DRY_RUN_TARGETS]
        sizes = sorted(set(self.lengths))
        # Probes that read probe_index are built for several probe numbers.
        indexes = range(8) if probe_index_nodes(self.protocol['probes']) else (0,)
        failures, empties, odd, lengths = {}, {}, {}, {}
        logger = logging.getLogger('uprobe.core.gen.probe')
        level = logger.level
        logger.setLevel(logging.CRITICAL)
        try:
            dag = DAG()
            dag.from_config({'probes': self.protocol['probes']})
            order, _ = _topological_order(dag.nodes)
            for size in sizes:
                region = _test_sequence(size)
                target_parts = {}
                if self.layout:
                    offset = 0
                    for name, part_size in zip(self.layout[0], self.layout[1][size][0]):
                        target_parts[name] = region[offset:offset + part_size]
                        offset += part_size
                for target, probe_index in ((t, i) for t in targets for i in indexes):
                    context = {'target_region': region, 'target': target, 'encoding': self.protocol['encoding'],
                               'target_parts': target_parts, PROBE_INDEX: probe_index}
                    for node in order:
                        node.done, node.result = False, None
                    for node in order:
                        if not all(dep.done for dep in node.deps):
                            continue
                        try:
                            node.build(context)
                        except Exception as exc:
                            failures.setdefault((node.name, f'{type(exc).__name__}: {exc}'), []).append((target, size))
                            continue
                        result = str(node.result)
                        self.built_lengths[node.name] = min(len(result), self.built_lengths.get(node.name, len(result)))
                        is_expr = self.nodes[node.name][0][-1] == 'expr'
                        if is_expr and not result and all(dep.result for dep in node.deps):
                            empties.setdefault(node.name, (target, size))  # root causes only
                        elif is_expr and set(result.upper()) - IUPAC:
                            odd.setdefault(node.name, result)
                        declared = self.nodes[node.name][1].get('length')
                        if type(declared) is int and len(result) != declared:
                            lengths.setdefault(node.name, (declared, len(result), size))
        except Exception as exc:
            self.warning(('probes',), f'probe construction could not be tested: {exc}',
                         'run `uprobe run` and check its log')
            return
        finally:
            logger.setLevel(level)
        total = len(targets) * len(sizes) * len(indexes)
        for (name, message), cases in failures.items():
            target, size = cases[0]
            self.error(self.nodes[name][0], f'fails when probes are built ({len(cases)}/{total} test cases, e.g. '
                       f'target {target!r} with a {size} nt target_region): {message}',
                       _failure_hint(message, self.lengths))
        for name, (target, size) in empties.items():
            self.warning(self.nodes[name][0], f'evaluates to an empty sequence (e.g. target {target!r} with a '
                         f'{size} nt target_region)', 'check the slice bounds against extracts.target_region.length')
        for name, result in odd.items():
            self.warning(self.nodes[name][0], f'evaluates to {result[:40]!r}, which is not a nucleotide sequence',
                         "return a sequence, e.g. encoding[target]['BC1'] rather than encoding[target]")
        for name, (declared, actual, size) in lengths.items():
            self.warning(self.nodes[name][0][:-1] + ('length',), f'declared length {declared} but the part builds '
                         f'{actual} nt (with a {size} nt target_region)', 'fix the expression or the declared length')

    def check_alignment_lengths(self, protocol, auto):
        """Warn when a bowtie2 mapped_genes attribute aligns sequences too short for its off-target search."""
        attributes = protocol.get('attributes')
        for name, attr in (attributes.items() if isinstance(attributes, dict) else []):
            if not isinstance(attr, dict) or str(attr.get('type', '')).lower() != 'mapped_genes' \
                    or attr.get('aligner') != 'bowtie2' or not isinstance(attr.get('target'), str):
                continue
            target = attr['target'].replace(':', '.')
            shortest = self.lengths[0] if target == 'target_region' and self.lengths else self.built_lengths.get(target)
            if not shortest or shortest > BOWTIE2_INCOMPLETE_MAX_LENGTH:
                continue
            path = ('attributes',) if auto else ('attributes', name, 'aligner')
            generated = ' (generated because `attributes` is empty)' if auto else ''
            if shortest < BOWTIE2_LOCAL_MIN_LENGTH:
                problem = (f'whose local mode cannot align sequences shorter than {BOWTIE2_LOCAL_MIN_LENGTH} nt: they '
                           f'count 0 genes, so filters on {name} keep them')
            else:
                problem = ('a seed heuristic that misses many off-target sites of sequences this short: it flagged '
                           'only 69-88% of 25-30-nt windows of human mRNAs that have a site within 2 substitutions '
                           'in another gene')
            if auto:
                fix = (f"install Bowtie 1 ({BOWTIE_INSTALL}) and add bowtie to the genome's align_index: the "
                       f'generated {name} then uses aligner: bowtie (every site within 2 substitutions); or define '
                       'attributes yourself')
            else:
                fix = ('aligner: bowtie   # every site within max_mismatches substitutions\n'
                       'max_mismatches: 2\n'
                       'strand: sense     # sense for target_region, antisense for a probe binding part\n'
                       f"and add bowtie to the genome's align_index (Bowtie 1: {BOWTIE_INSTALL})")
            self.warning(path, f'{name}{generated} aligns {target} ({shortest} nt) with bowtie2, {problem}', fix)


def validate_protocol(protocol, genomes=None, check_targets: bool = False, *,
                      files: T.Optional[dict] = None, lines: T.Optional[dict] = None) -> T.List[Finding]:
    """Check a loaded protocol (and optional genome config); return findings, errors first.

    ``files`` and ``lines`` map 'protocol'/'genomes' to the file name and to a
    {path tuple: line} index; validate_files fills them in.
    """
    findings = _Validator(protocol, genomes, check_targets, files, lines).run()
    return sorted(findings, key=lambda f: f.severity != 'error')


def _line_index(text):
    """Map key paths to 1-based lines; also return duplicate keys as (path, line, first line)."""
    index, duplicates = {}, []
    try:
        root = yaml.compose(text, Loader=yaml.SafeLoader)
    except yaml.YAMLError:
        return index, duplicates

    def walk(node, path):
        if isinstance(node, yaml.MappingNode):
            seen = {}
            for key_node, value_node in node.value:
                if not isinstance(key_node, yaml.ScalarNode):
                    continue
                key, line = key_node.value, key_node.start_mark.line + 1
                child = path + (key,)
                if key in seen:
                    duplicates.append((child, line, seen[key]))
                seen[key] = line
                index[child] = line
                walk(value_node, child)
        elif isinstance(node, yaml.SequenceNode):
            for i, item in enumerate(node.value):
                index[path + (i,)] = item.start_mark.line + 1
                walk(item, path + (i,))

    walk(root, ())
    return index, duplicates


def _load_yaml(path: Path, doc: str):
    """Return (data, line index, findings) for one YAML file."""
    label = str(path)
    try:
        text = path.read_text(encoding='utf-8')
    except (OSError, UnicodeDecodeError) as exc:
        return None, {}, [Finding('error', '', f'cannot read the file: {exc}', 'pass a UTF-8 YAML file',
                                  label, None, (), doc)]
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        mark = getattr(exc, 'problem_mark', None)
        line = mark.line + 1 if mark else None
        problem = getattr(exc, 'problem', None) or str(exc)
        context = getattr(exc, 'context', None)
        source_line = text.splitlines()[mark.line] if mark and mark.line < len(text.splitlines()) else ''
        hint = ('indent with spaces, not tabs' if '\t' in source_line else
                'check the indentation and quote values containing ": ", " #", "{" or "[", '
                'e.g. expr: "encoding[target][\'BC1\']"')
        return None, {}, [Finding('error', '', f'YAML syntax error: {problem}' + (f' ({context})' if context else ''),
                                  hint, label, line, (), doc)]
    if data is None:
        return None, {}, [Finding('error', '', 'the file is empty', 'see README.md for an example',
                                  label, None, (), doc)]
    index, duplicates = _line_index(text)
    findings = [Finding('warning', format_path(key_path), f'duplicate key {key_path[-1]!r} (first on line {first}); '
                        'YAML keeps only the last value', 'merge the two entries or rename one',
                        label, line, key_path, doc)
                for key_path, line, first in duplicates]
    return data, index, findings


def validate_files(protocol_path, genomes_path=None, check_targets: bool = False) -> T.List[Finding]:
    """Load and check a protocol YAML (and optional genomes YAML)."""
    files, lines, findings = {}, {}, []
    protocol, index, problems = _load_yaml(Path(protocol_path), 'protocol')
    files['protocol'], lines['protocol'] = str(protocol_path), index
    findings += problems
    genomes = None
    if genomes_path is not None:
        genomes, index, problems = _load_yaml(Path(genomes_path), 'genomes')
        files['genomes'], lines['genomes'] = str(genomes_path), index
        findings += problems
    if protocol is None:
        return sorted(findings, key=lambda f: f.severity != 'error')
    loaded_genomes = genomes if genomes_path is not None and not any(
        f.doc == 'genomes' and f.severity == 'error' for f in findings) else None
    findings += validate_protocol(protocol, loaded_genomes, check_targets, files=files, lines=lines)
    return sorted(findings, key=lambda f: f.severity != 'error')


def summarize(findings: T.List[Finding]) -> dict:
    errors = sum(f.severity == 'error' for f in findings)
    return {'valid': errors == 0, 'errors': errors, 'warnings': len(findings) - errors}


def render_json(findings: T.List[Finding]) -> str:
    return json.dumps({'summary': summarize(findings), 'findings': [f.to_dict() for f in findings]},
                      indent=2, ensure_ascii=False)


def render_text(findings: T.List[Finding]) -> str:
    out = []
    for f in findings:
        where = (f.file or f.doc) + (f':{f.line}' if f.line else '')
        out.append(f'{where}: {f.severity}: {f.path or "<root>"}: {f.message}')
        if f.suggestion:
            first, *rest = f.suggestion.splitlines()
            out.append(f'  fix: {first}')
            out.extend(f'       {line}' for line in rest)
    s = summarize(findings)
    out.append(f"{'OK' if s['valid'] else 'FAILED'}: {s['errors']} error(s), {s['warnings']} warning(s)")
    return '\n'.join(out)
