"""Where the indexes of a genomes.yaml entry live.

By default next to the entry's FASTA: <aligner>_genome/<stem>, <aligner>_transcript/<stem>,
jf_genome/ and transcript.fa. Two optional keys of the entry move them::

    index_dir: /fast/disk/indexes/hg38      # the folders above go here instead of next to the FASTA
    indexes:                                # an index that already exists: its prefix, as given to bowtie
      bowtie_transcript: /data/gencode/bowtie/transcripts

``indexes`` keys are <aligner>_<kind>, with kind genome or transcript. An index listed there is
used as it is and never built; a listed index that does not exist is an error.
"""
from pathlib import Path

KINDS = ('genome', 'transcript')


def index_root(genome: dict) -> Path:
    """The directory holding the derived folders and transcript.fa of the entry."""
    if genome.get('index_dir'):
        return Path(str(genome['index_dir']))
    return Path(str(genome['fasta'])).parent


def pointed_index(genome: dict, aligner: str, kind: str):
    """The prefix of an existing index listed under ``indexes``, or None."""
    value = (genome.get('indexes') or {}).get(f'{aligner}_{kind}')
    return Path(str(value)) if value else None


def index_prefix(genome: dict, aligner: str, kind: str) -> Path:
    """The prefix of the ``kind`` index of ``aligner`` (genome or transcript)."""
    return pointed_index(genome, aligner, kind) or \
        index_root(genome) / f'{aligner}_{kind}' / Path(str(genome['fasta'])).stem
