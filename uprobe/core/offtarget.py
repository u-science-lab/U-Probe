"""Settings of ``mapped_genes``, the transcriptome off-target search.

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
"""

MAPPED_GENES_ALIGNERS = ('bowtie2', 'bowtie')
# Command-line options of both Bowtie and Bowtie2 for each strand setting.
STRAND_OPTIONS = {'both': [], 'sense': ['--norc'], 'antisense': ['--nofw']}
STRANDS = tuple(STRAND_OPTIONS)
DEFAULT_MAX_MISMATCHES = 2
MAX_MISMATCHES_LIMIT = 3  # bowtie -v accepts 0-3
# Shortest sequence Bowtie2's local mode (default --score-min G,20,8) can align.
BOWTIE2_LOCAL_MIN_LENGTH = 22


def short_for_bowtie2(attr_name: str, short: int, total: int) -> str:
    """Message for sequences that a bowtie2 mapped_genes attribute cannot align."""
    return (f'{attr_name}: {short} of {total} sequences are shorter than {BOWTIE2_LOCAL_MIN_LENGTH} nt, which '
            f'Bowtie2\'s local mode cannot align (its minimum score exceeds a perfect match), so they count 0 '
            f'genes; use aligner: bowtie for an exhaustive search of short sequences')


def mapped_genes_settings(attr: dict, where: str = '') -> dict:
    """Return {aligner, strand, max_mismatches} of a mapped_genes attribute.

    Raises ValueError starting with the option name (after ``where.``);
    ``max_mismatches`` is None for bowtie2.
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
    return {'aligner': aligner, 'strand': strand, 'max_mismatches': max_mismatches}
