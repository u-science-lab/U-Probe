import sys
import os
import contextlib
from pathlib import Path
import logging
import click
import yaml
import copy
import pandas as pd


from .api import UProbeAPI
from .utils import get_logger
from uprobe import __version__


log = get_logger(__name__)

# Exit status of `uprobe run`; click exits with 2 on usage errors.
EXIT_ERROR = 1
EXIT_NO_PROBES = 3


def _extract_probe_targets(probes_config: dict) -> list:
    """
    Extract all probe and part targets from probes configuration.
    Returns list of (target_name, is_probe_level) tuples.
    Examples: [('probe_1', True), ('probe_1.part1', False), ('probe_1.part2', False)]
    """
    if not isinstance(probes_config, dict):
        raise ValueError(
            "Invalid protocol: `probes` must be a YAML mapping (dict) like "
            "`probes: { probe_1: {template: ..., parts: {...}} }`. "
            "Do NOT use list-style probes."
        )

    targets = []
    
    def traverse_parts(prefix: str, config: dict, is_top_level: bool = False):
        """Recursively traverse probe parts."""
        if is_top_level:
            targets.append((prefix, True))  # probe level
        
        if 'parts' in config:
            if not isinstance(config['parts'], dict):
                raise ValueError(
                    f"Invalid protocol: `{prefix}.parts` must be a mapping (dict). "
                    "Do NOT use list-style parts; use `parts: {part1: {...}, part2: {...}}`."
                )
            for part_name, part_config in config['parts'].items():
                part_full_name = f"{prefix}.{part_name}"
                targets.append((part_full_name, False))  # part level
                # Recursively handle nested parts
                if 'parts' in part_config:
                    traverse_parts(part_full_name, part_config, is_top_level=False)
    
    for probe_name, probe_config in probes_config.items():
        traverse_parts(probe_name, probe_config, is_top_level=True)
    
    return targets


def _default_mapped_genes(protocol_config: dict, genome) -> dict:
    """The generated target_mapped_genes: Bowtie 1 when the genome can use it, else Bowtie2.

    Bowtie2's local mode aligns no sequence shorter than 22 nt and misses
    sites of somewhat longer ones; Bowtie 1 finds every site within 2
    substitutions. target_region is the RNA sequence, which the probes bind
    where it occurs: strand sense.
    """
    from .offtarget import BOWTIE2_LOCAL_MIN_LENGTH, DEFAULT_MAX_MISMATCHES, bowtie_unavailable
    from .sampling import normalize_sampling
    reason = bowtie_unavailable(genome)
    if not reason:
        log.info(f"Auto-generated target_mapped_genes: aligner bowtie, every transcript site within "
                 f"{DEFAULT_MAX_MISMATCHES} substitutions (strand sense)")
        return {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                'max_mismatches': DEFAULT_MAX_MISMATCHES, 'strand': 'sense'}
    message = (f"Auto-generated target_mapped_genes: aligner bowtie2, because {reason}. Install Bowtie 1 "
               "(conda install -c bioconda bowtie) and add bowtie to the genome's align_index for an exhaustive "
               "off-target search")
    try:
        shortest = normalize_sampling(protocol_config['extracts']['target_region'].get('length'))[0]
    except (TypeError, ValueError):
        shortest = None
    if shortest is not None and shortest < BOWTIE2_LOCAL_MIN_LENGTH:
        log.warning(f"{message}: Bowtie2 aligns no sequence shorter than {BOWTIE2_LOCAL_MIN_LENGTH} nt, so the "
                    f"{shortest}-nt windows count no off-target gene.")
    else:
        log.info(message + ".")
    return {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie2'}


# The generated kmer_count attributes of DNA mode (Jellyfish, 35-mers).
DNA_KMER_COUNT = {'type': 'kmer_count', 'kmer_len': 35, 'size': '1G', 'aligner': 'jellyfish'}


def _generate_default_attributes(protocol_config: dict, genome: dict = None) -> dict:
    """
    Generate default attributes based on probes structure and mode (DNA/RNA).
    DNA mode focuses on probe parts; RNA mode includes target_region + all probes/parts.
    `genome` is the genome entry of the run, if known: it decides the aligner of target_mapped_genes.
    """
    source = protocol_config['extracts']['target_region'].get('source', 'genome')
    is_dna_mode = (source == 'genome')
    probes_config = protocol_config.get('probes', {})
    
    attributes = {}
    
    # DNA mode: minimal or no target_region attributes (focus on probe parts)
    # RNA mode: comprehensive target_region attributes
    if not is_dna_mode:
        # RNA mode: add full target_region attributes
        attributes['target_gc'] = {
            'target': 'target_region',
            'type': 'gc_content'
        }
        attributes['target_tm'] = {
            'target': 'target_region',
            'type': 'annealing_temperature'
        }
        attributes['target_fold'] = {
            'target': 'target_region',
            'type': 'fold_score'
        }
        attributes['target_self_match'] = {
            'target': 'target_region',
            'type': 'self_match'
        }
        attributes['target_mapped_genes'] = _default_mapped_genes(protocol_config, genome)
    else:
        # The genomic window is the sequence that binds the genome, so its k-mer
        # count is the repeat measure the generated filter target_kmer_count
        # <= 100 reads. Jellyfish counts canonical k-mers, so it equals that of
        # a binding part rc(target_region); constant parts (primers, barcodes)
        # do not bind the genome and must not decide which windows pass.
        attributes['target_kmer_count'] = dict(DNA_KMER_COUNT, target='target_region')

    # Extract probe targets and generate attributes
    probe_targets = _extract_probe_targets(probes_config)
    
    for target_name, is_probe_level in probe_targets:
        # Generate safe attribute name (replace dots with underscores)
        safe_name = target_name.replace('.', '_')
        
        # DNA mode: focus on probe parts only
        # RNA mode: include both probe-level and part-level
        if is_dna_mode and is_probe_level:
            # DNA: skip probe-level attributes, only process parts
            continue
        
        # Common attributes for probes/parts
        attributes[f'{safe_name}_gc'] = {
            'target': target_name,
            'type': 'gc_content'
        }
        attributes[f'{safe_name}_tm'] = {
            'target': target_name,
            'type': 'annealing_temperature'
        }
        attributes[f'{safe_name}_fold'] = {
            'target': target_name,
            'type': 'fold_score'
        }
        attributes[f'{safe_name}_self_match'] = {
            'target': target_name,
            'type': 'self_match'
        }
        
        # DNA mode: add specificity attributes to probe parts
        if is_dna_mode and not is_probe_level:
            # Add mapped_sites and kmer_count for DNA probe parts
            attributes[f'{safe_name}_mapped_sites'] = {
                'target': target_name,
                'type': 'mapped_sites',
                'aligner': 'bowtie2'
            }
            attributes[f'{safe_name}_kmer_count'] = dict(DNA_KMER_COUNT, target=target_name)
    
    return attributes


def _generate_default_summary(protocol_config: dict, attributes: dict) -> dict:
    """
    Generate default summary configuration based on mode (DNA/RNA) and attributes.
    DNA: report_name=dna_report, includes probe part attributes
    RNA: report_name=rna_report, includes target + probe-level attributes
    """
    source = protocol_config['extracts']['target_region'].get('source', 'genome')
    is_dna_mode = (source == 'genome')
    
    summary = {
        'report_name': 'dna_report' if is_dna_mode else 'rna_report',
        'attributes': []
    }
    
    # Select attributes for summary based on mode
    for attr_name in attributes.keys():
        if is_dna_mode:
            # DNA: the window's k-mer count and part-level attributes (gc, tm, kmer_count)
            if attr_name == 'target_kmer_count' or (
                    '_part' in attr_name and any(attr_name.endswith(x) for x in ['_gc', '_tm', '_kmer_count'])):
                summary['attributes'].append(attr_name)
        else:
            # RNA: include target + probe-level attributes (not parts)
            if attr_name.startswith('target_'):
                # Include target attributes
                if any(attr_name.endswith(x) for x in ['_gc', '_tm', '_fold', '_self_match', '_mapped_genes']):
                    summary['attributes'].append(attr_name)
            elif '_part' not in attr_name:
                # Include probe-level attributes (not parts)
                if any(attr_name.endswith(x) for x in ['_gc', '_tm', '_fold']):
                    summary['attributes'].append(attr_name)
            else:
                # Include key part attributes for RNA as well
                if any(attr_name.endswith(x) for x in ['_tm']):
                    summary['attributes'].append(attr_name)
    
    return summary


def _number(value: float) -> str:
    """A number as filter conditions read it: plain decimal notation."""
    return f'{value:f}'.rstrip('0').rstrip('.')


def _tm_condition(protocol_config: dict, attr_name: str) -> str:
    """Tm filter: 50-90 °C, or [T, T + 40] °C at the declared hybridization temperature T."""
    from .hybridization import default_tm_window
    low, high, _ = default_tm_window(protocol_config, attr_name)
    return f'{attr_name} >= {_number(low)} & {attr_name} <= {_number(high)}'


def _log_tm_windows(protocol_config: dict) -> None:
    """Log the window of each auto-generated Tm filter and why it was chosen."""
    from .hybridization import default_tm_window
    reasons = {'protocol': 'hybridization.temperature_c to +40 °C',
               'attribute': "the attribute's hybridization.temperature_c to +40 °C",
               'default': 'default 37 °C to +40 °C; hybridization.temperature_c is not set',
               None: 'legacy window for the legacy Tm conditions'}
    windows = {}
    for name in protocol_config['post_process'].get('filters', {}):
        if name.endswith('_tm'):
            windows.setdefault(default_tm_window(protocol_config, name), []).append(name)
    for (low, high, source), names in windows.items():
        log.info(f"Auto-generated Tm filter {_number(low)}-{_number(high)} °C ({reasons[source]}): {', '.join(names)}")
        if source == 'default':
            log.warning("hybridization sets Tm conditions but no temperature_c, so the auto-generated Tm "
                        f"filters assume 37 °C ({_number(low)}-{_number(high)} °C); declare your "
                        "hybridization temperature as hybridization.temperature_c.")


# Generated score attributes ranked lowest first: the least structured, least
# self-complementary and least repetitive candidates come first.
LOWER_IS_BETTER_SUFFIXES = ('_fold', '_self_match', '_kmer_count')


def _generate_default_post_process(protocol_config: dict, attributes: dict) -> dict:
    """
    Generate default post_process configuration based on attributes.
    Tm filters are 50-90 °C, or [T, T + 40] °C at the declared hybridization temperature T.
    """
    source = protocol_config['extracts']['target_region'].get('source', 'genome')
    is_dna_mode = (source == 'genome')
    tm_protocol = {**protocol_config, 'attributes': attributes}
    
    post_process = {
        'filters': {},
        'sorts': {
            'is_ascending': [],
            'is_descending': []
        }
    }
    
    # Generate filters for target_region attributes
    if 'target_gc' in attributes:
        post_process['filters']['target_gc'] = {
            'condition': 'target_gc >= 0.2 & target_gc <= 0.8'  # GC is 0-1 fraction
        }
    
    if 'target_tm' in attributes:
        post_process['filters']['target_tm'] = {
            'condition': _tm_condition(tm_protocol, 'target_tm')  # Tm in Celsius
        }
    
    # Mode-specific filters: genome repeats (DNA), off-target genes (RNA)
    if is_dna_mode:
        if 'target_kmer_count' in attributes:
            post_process['filters']['target_kmer_count'] = {
                'condition': 'target_kmer_count <= 100'
            }
    else:
        if 'target_mapped_genes' in attributes:
            post_process['filters']['target_mapped_genes'] = {
                'condition': 'target_mapped_genes <= 4'
            }
    
    # Generate filters for probe/part attributes (only Tm filters for parts)
    for attr_name in attributes.keys():
        if attr_name.startswith('target_'):
            continue
        
        # Add Tm filters for all probes/parts
        if attr_name.endswith('_tm'):
            post_process['filters'][attr_name] = {
                'condition': _tm_condition(tm_protocol, attr_name)
            }
        
        # Add GC filters for probe-level only (not parts)
        if attr_name.endswith('_gc') and '.' not in attributes[attr_name]['target']:
            post_process['filters'][attr_name] = {
                'condition': f'{attr_name} >= 0.2 & {attr_name} <= 0.8'
            }
    
    # Generate sorts, all ascending: GC and Tm first, then the scores where
    # lower is better: fewer mapped genes, less structure (fold_score is -MFE,
    # >= 0), less self-complementarity, lower genome k-mer counts.
    for attr_name in attributes.keys():
        if (attr_name.endswith('_gc') or attr_name.endswith('_tm')
                or attr_name.endswith('_mapped_genes')):
            post_process['sorts']['is_ascending'].append(attr_name)
    for attr_name in attributes.keys():
        if any(attr_name.endswith(suffix) for suffix in LOWER_IS_BETTER_SUFFIXES):
            post_process['sorts']['is_ascending'].append(attr_name)

    return post_process


def _is_blank_step(value) -> bool:
    """True for a post_process step left as a placeholder: empty, or every option blank."""
    if isinstance(value, dict) and value:
        return all(option is None or option == '' for option in value.values())
    return value in ({}, [], None)


def _validate_and_normalize_protocol(protocol_config: dict, genome: dict = None) -> dict:
    """
    Validate protocol configuration and inject defaults for optional fields.
    Dynamically generates attributes and post_process based on probes structure.
    `genome`, the genome entry of the run (genomes.yaml), decides whether a
    generated off-target attribute can use Bowtie 1; without it, Bowtie2.
    """
    # 1. Check required top-level keys
    required_keys = ['genome', 'targets', 'extracts', 'encoding', 'probes']
    missing_keys = [k for k in required_keys if k not in protocol_config]
    if missing_keys:
        raise ValueError(f"Protocol is missing required keys: {', '.join(missing_keys)}")
    
    # Check extracts.target_region
    if 'target_region' not in protocol_config.get('extracts', {}):
        raise ValueError("Protocol must define 'extracts.target_region'")
    
    from .sampling import normalize_sampling
    region = protocol_config['extracts']['target_region']
    _, _, region['step'] = normalize_sampling(region.get('length'), region.get('step'), region.get('overlap'))
    region.pop('overlap', None)
    if 'transcripts' in region:
        from .gen.fun import transcripts_setting
        try:
            region['transcripts'] = transcripts_setting(region['transcripts'])
        except ValueError as exc:
            raise ValueError(f'extracts.target_region.{exc}') from None
        if region.get('source', 'genome') == 'genome' and region['transcripts'] != 'all':
            log.warning('extracts.target_region.transcripts applies to gene targets (source mRNA, exon, CDS or '
                        'UTR); it is ignored for source: genome')
    if 'layout' in region:
        from .layout import compile_layout
        compile_layout(region['layout'], region['length'])
        if 'target_parts' in protocol_config.get('probes', {}):
            raise ValueError('target_parts is reserved for layout references')

    from .layout import validate_layout_references
    validate_layout_references(protocol_config)

    # 2. Determine mode (DNA vs RNA)
    source = protocol_config['extracts']['target_region'].get('source', 'genome')
    is_dna_mode = (source == 'genome')

    # 2.1 Clean up user-provided post_process (remove empty stubs and DNA-only keys in RNA mode)
    # Agents/users sometimes write placeholders like:
    #   post_process:
    #     filters: {}
    #     avoid_otp: {}
    #     equal_space: {}
    # These should be treated as "missing" so defaults can be generated, and RNA mode
    # must not carry the DNA-only avoid_otp step (equal_space works in both modes).
    pp = protocol_config.get('post_process')
    if isinstance(pp, dict):
        if not is_dna_mode:
            pp.pop('avoid_otp', None)
        # Drop empty/None/list stubs at top-level, and steps whose options are
        # all blank, e.g. `remove_overlap: {location_interval: }`
        pp = {k: v for k, v in pp.items() if not _is_blank_step(v)}
        protocol_config['post_process'] = pp
    
    # 3. Auto-complete attributes if missing (dynamic generation based on probes structure)
    if 'attributes' not in protocol_config or not protocol_config['attributes']:
        log.info(f"Attributes missing or empty. Auto-generating from probes structure for mode: {'DNA' if is_dna_mode else 'RNA'}")
        protocol_config['attributes'] = _generate_default_attributes(protocol_config, genome)
    
    # 3.1 Validate hybridization conditions (protocol block and per-attribute overrides)
    from .hybridization import validate_protocol_hybridization
    validate_protocol_hybridization(protocol_config)

    # 3.2 Validate the gene_expression block; its table is read when mapped_genes attributes run.
    if protocol_config.get('gene_expression') is not None:
        from .gene_expression import validate_gene_expression
        validate_gene_expression(protocol_config['gene_expression'])
    
    # 4. Auto-complete post_process if missing (dynamic generation based on attributes)
    if 'post_process' not in protocol_config or not protocol_config['post_process']:
        log.info("Post-process configuration missing or empty. Auto-generating from attributes.")
        protocol_config['post_process'] = _generate_default_post_process(
            protocol_config, 
            protocol_config['attributes']
        )
        _log_tm_windows(protocol_config)
    
    # 4.1 Check off_target_pool before the slow steps: it reads mapped_genes attributes.
    pool = (protocol_config.get('post_process') or {}).get('off_target_pool')
    if pool:
        from .process.off_target_pool import check_pool_attributes, pool_settings
        check_pool_attributes(pool_settings(pool)[0], protocol_config['attributes'])

    # 5. Auto-complete summary if missing (dynamic generation based on mode and attributes)
    if 'summary' not in protocol_config or not protocol_config['summary']:
        log.info(f"Summary configuration missing or empty. Auto-generating with report_name: {'dna_report' if is_dna_mode else 'rna_report'}")
        protocol_config['summary'] = _generate_default_summary(
            protocol_config,
            protocol_config['attributes']
        )
    
    return protocol_config


@click.group()
@click.version_option(
    __version__,
    '--version',
    '-V',
    prog_name='uprobe',
    message='U-Probe %(version)s',
)
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging.')
@click.option('--quiet', '-q', is_flag=True, help='Suppress all output except errors.')
@click.pass_context
def cli(ctx, verbose, quiet):
    """
    U-Probe: Universal Probe Design Tool
    
    A powerful and flexible Python-based tool for designing custom DNA or RNA probes
    for various molecular biology applications.
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    elif verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)
    ctx.ensure_object(dict)


@cli.command()
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True), 
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--output', '-o', default='./results', type=click.Path(),
              help='Output directory. [default: ./results]')
@click.option('--raw', is_flag=True,
              help='Save unfiltered raw probe data.')
@click.option('--continue-invalid', is_flag=True,
              help='Skip targets that are not in the genome annotation and design the others. '
                   'Without it, any such target stops the run before designing.')
@click.option('--threads', '-t', default=10, type=click.IntRange(min=1),
              help='Threads for index building and for the aligner and Jellyfish steps of every attribute '
                   'that sets no threads of its own. [default: 10]')
@click.option('--json', 'json_output', is_flag=True,
              help='Print the run summary (summary.json) to stdout; logs and tool output go to stderr.')
def run(protocol, genomes, output, raw, continue_invalid, threads, json_output):
    """
    Execute the complete probe design workflow.
    
    This command runs the entire pipeline from genome index building to final probe generation.
    Every run writes <output>/summary.json and the final table as <name>_<time>_probes.tsv.
    A gene target missing from the GTF stops the run before designing (status error),
    unless --continue-invalid skips it.

    \b
    Exit status (with 0, 1 and 3, summary.json and the --json output are written first):
      0  probes were designed (status ok or partial)
      1  the run failed (status error), e.g. a target not in the GTF without --continue-invalid
      2  invalid command-line usage (missing or bad options)
      3  the run finished without any probe (status no_probes,
         no_probes_constructed, no_target_sequences or no_valid_targets)
    """
    uprobe = None
    protocol_config = None
    error = None
    # --json: keep stdout clean for the summary; logs and tool output go to stderr.
    stdout_guard = contextlib.ExitStack()
    if json_output:
        stdout_guard.enter_context(_stdout_to_stderr())
    try:
        log.info("Starting U-Probe complete workflow...")
        
        # Load and validate protocol
        protocol_path = Path(protocol)
        with open(protocol_path, 'r', encoding='utf-8') as f:
            protocol_config = yaml.safe_load(f) or {}
            
        if not isinstance(protocol_config, dict):
            raise ValueError("Protocol YAML must be a mapping (dict)")
            
        # Normalize protocol (inject defaults); the genome entry decides the
        # aligner of a generated off-target attribute.
        protocol_config = _validate_and_normalize_protocol(protocol_config, _genome_entry(genomes, protocol_config))
        
        uprobe = UProbeAPI(
            protocol_config=protocol_config,
            genomes_config=Path(genomes),
            output_dir=Path(output)
        )
        result_df = uprobe.run_workflow(
            raw_csv=raw,
            continue_on_invalid_targets=continue_invalid,
            threads=threads
        )
        if not result_df.empty:
            log.info(f"Workflow completed successfully with {len(result_df)} final probes!")
    except Exception as e:
        log.error(f"Workflow failed: {e}")
        error = e
    finally:
        stdout_guard.close()
    summary = getattr(uprobe, 'summary', None)
    if summary is None:
        summary = _write_failure_summary(protocol_config, error, Path(output))
    _log_target_counts(summary)
    if json_output:
        from .report.run_summary import summary_json
        click.echo(summary_json(summary))
    if error is not None:
        sys.exit(EXIT_ERROR)
    from .report.run_summary import NO_PROBE_STATUSES
    if summary.get('status') in NO_PROBE_STATUSES:
        log.error(f"No probes designed (status {summary['status']}): {summary.get('message')} "
                  f"See {Path(output) / 'summary.json'}.")
        sys.exit(EXIT_NO_PROBES)


def _genome_entry(genomes_path, protocol_config: dict):
    """The protocol's genome entry in the genome config, or None (UProbeAPI reports a missing one)."""
    with open(genomes_path, 'r', encoding='utf-8') as f:
        genomes = yaml.safe_load(f)
    name = protocol_config.get('genome')
    entry = genomes.get(name) if isinstance(genomes, dict) and isinstance(name, str) else None
    return entry if isinstance(entry, dict) else None


@contextlib.contextmanager
def _stdout_to_stderr():
    """Send everything written to stdout, including tool subprocesses, to stderr."""
    try:
        stdout_fd, stderr_fd = sys.stdout.fileno(), sys.stderr.fileno()
    except (AttributeError, OSError, ValueError):
        # No real descriptors (e.g. click's CliRunner): redirect Python writes only.
        with contextlib.redirect_stdout(sys.stderr):
            yield
        return
    sys.stdout.flush()
    saved_fd = os.dup(stdout_fd)
    os.dup2(stderr_fd, stdout_fd)
    try:
        with contextlib.redirect_stdout(sys.stderr):
            yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, stdout_fd)
        os.close(saved_fd)


def _log_target_counts(summary: dict) -> None:
    """Final probes per target, with the number requested by equal_space and any shortfall."""
    targets = summary.get('targets') or {}
    if summary.get('status') == 'error' or not targets:
        return
    counts = ', '.join(f"{name} {entry['final']}" + (f"/{entry['requested_probes']}" if entry.get('requested_probes') else '')
                       for name, entry in targets.items())
    requested = any(entry.get('requested_probes') for entry in targets.values())
    log.info(f"Final probes per target{' (final/requested)' if requested else ''}: {counts}")
    for warning in summary.get('warnings') or []:
        if warning.get('code') == 'fewer_probes_than_requested':
            short = ', '.join(f"{name} {c['final']} of {c['requested']}" for name, c in warning['counts'].items())
            log.warning(f"{warning['message']} {short}")


def _write_failure_summary(protocol_config, error, output_dir: Path) -> dict:
    """Summarize a run that failed before the workflow could write summary.json."""
    from .report.run_summary import SUMMARY_FILENAME, build_run_summary, write_run_summary
    summary = build_run_summary(protocol_config or {}, {'step': 'load_inputs',
                                                        'error': error or RuntimeError('Run summary unavailable')})
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_run_summary(summary, output_dir / SUMMARY_FILENAME)
    except OSError as e:
        log.error(f"Could not write run summary: {e}")
    return summary


@cli.command()
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True, dir_okay=False),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', type=click.Path(exists=True, dir_okay=False),
              help='Path to genome configuration file (YAML). Enables genome and aligner checks.')
@click.option('--json', 'as_json', is_flag=True,
              help='Print the findings and a summary as JSON.')
@click.option('--check-targets', is_flag=True,
              help='Also look up genome targets in the GTF (gene names) or FASTA (regions). Requires --genomes.')
def validate(protocol, genomes, as_json, check_targets):
    """
    Check a protocol before running it.

    Reports every problem with its location in the YAML and a suggested fix.
    Exits with status 1 when there are errors; warnings alone exit with 0.
    """
    from .validate import validate_files, render_json, render_text, summarize
    if check_targets and not genomes:
        raise click.UsageError("'--check-targets' requires '--genomes'.")
    findings = validate_files(Path(protocol), Path(genomes) if genomes else None, check_targets=check_targets)
    click.echo(render_json(findings) if as_json else render_text(findings))
    sys.exit(0 if summarize(findings)['valid'] else 1)


@cli.command(name='build-index')
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--threads', '-t', default=10, type=int,
              help='Number of threads for index building. [default: 10]')
def build_index(protocol, genomes, threads):
    """
    Build genome index for alignment tools.
    
    Creates necessary indices for tools like Bowtie2 and BLAST as specified
    in the genome configuration.
    """
    try:
        log.info("Building genome index...")
        uprobe = UProbeAPI(
            protocol_config=Path(protocol),
            genomes_config=Path(genomes),
            output_dir=Path('./temp') 
        )
        
        uprobe.build_genome_index(threads=threads)
        log.info("Genome index building completed successfully!")
    except Exception as e:
        log.error(f"Index building failed: {e}")
        sys.exit(1)


@cli.command(name='validate-targets')
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--continue-invalid', is_flag=True,
              help='Exit with 0 when some targets are invalid but others are valid.')
def validate_targets(protocol, genomes, continue_invalid):
    """
    Validate target genes against the genome annotation.

    Checks if all target genes specified in the protocol exist in the GTF file.
    """
    try:
        log.info("Validating target genes...")
        uprobe = UProbeAPI(
            protocol_config=Path(protocol),
            genomes_config=Path(genomes),
            output_dir=Path('./temp')
        )
        valid = uprobe.validate_targets(continue_on_invalid=continue_invalid)
        if not valid:
            log.error("Target validation failed: no valid targets.")
            sys.exit(1)
        if uprobe.target_validation['invalid']:
            log.warning("Some targets are invalid; the others are valid (--continue-invalid).")
        else:
            log.info("All target genes are valid!")
    except Exception as e:
        log.error(f"Target validation failed: {e}")
        sys.exit(1)


@cli.command(name='generate-targets')
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--output', '-o', default='./results', type=click.Path(),
              help='Output directory. [default: ./results]')
@click.option('--continue-invalid', is_flag=True,
              help='Skip targets that are not in the genome annotation; without it, any such target stops.')
def generate_targets(protocol, genomes, output, continue_invalid):
    """
    Generate target region sequences from the genome.
    
    Extracts target sequences based on the extraction parameters in the protocol.
    """
    try:
        log.info("Generating target sequences...")
        uprobe = UProbeAPI(
            protocol_config=Path(protocol),
            genomes_config=Path(genomes),
            output_dir=Path(output)
        )
        if not uprobe.validate_targets(continue_on_invalid=continue_invalid):
            log.error("Target validation failed!")
            sys.exit(1)
        
        df_targets = uprobe.generate_target_seqs()
        
        if df_targets.empty:
            log.error("No target sequences generated!")
            sys.exit(1)
        else:
            log.info(f"Generated {len(df_targets)} target sequences successfully!")
            targets_file = Path(output) / "target_sequences.csv"
            Path(output).mkdir(parents=True, exist_ok=True)
            df_targets.to_csv(targets_file, index=False)
            log.info(f"Target sequences saved to {targets_file}")           
    except Exception as e:
        log.error(f"Target generation failed: {e}")
        sys.exit(1)


@cli.command(name='construct-probes')
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--targets', required=True, type=click.Path(exists=True),
              help='Path to target sequences CSV file.')
@click.option('--output', '-o', default='./results', type=click.Path(),
              help='Output directory. [default: ./results]')
def construct_probes(protocol, genomes, targets, output):
    """
    Construct probes from target sequences.
    
    Takes target sequences and constructs probes according to the probe design
    specifications in the protocol.
    """
    try:
        import pandas as pd
        
        log.info("Constructing probes...")
        uprobe = UProbeAPI(
            protocol_config=Path(protocol),
            genomes_config=Path(genomes),
            output_dir=Path(output)
        )
        df_targets = pd.read_csv(targets)        
        df_probes = uprobe.construct_probes(df_targets)       
        if df_probes.empty:
            log.error("No probes constructed!")
            sys.exit(1)
        else:
            log.info(f"Constructed {len(df_probes)} probes successfully!")
            out_dir = Path(output)
            out_dir.mkdir(parents=True, exist_ok=True)
            probes_file = out_dir / "constructed_probes.csv"
            df_probes.to_csv(probes_file, index=False)
            log.info(f"Constructed probes saved to {probes_file}")

            if len(df_targets) == len(df_probes):
                combined_file = out_dir / "constructed_probes_combined.csv"
                df_combined = pd.concat(
                    [df_targets.reset_index(drop=True), df_probes.reset_index(drop=True)],
                    axis=1,
                )
                df_combined.to_csv(combined_file, index=False)
                log.info(f"Combined target/probe table saved to {combined_file}")
            else:
                log.warning(
                    "Target/probe row counts differ; combined post-process input was not saved."
                )
    except Exception as e:
        log.error(f"Probe construction failed: {e}")
        sys.exit(1)


@cli.command(name='post-process')
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--probes', required=True, type=click.Path(exists=True),
              help='Path to probe data CSV file (combined targets and probes).')
@click.option('--output', '-o', default='./results', type=click.Path(),
              help='Output directory. [default: ./results]')
@click.option('--raw', is_flag=True,
              help='Save unfiltered raw probe data.')
@click.option('--threads', '-t', default=10, type=click.IntRange(min=1),
              help='Threads for the aligner and Jellyfish steps of every attribute that sets no threads '
                   'of its own. [default: 10]')
def post_process(protocol, genomes, probes, output, raw, threads):
    """
    Post-process probes (add attributes and apply filters).
    
    Adds various attributes to probes and applies filtering criteria
    as specified in the protocol.
    """
    try:
        log.info("Post-processing probes...")
        uprobe = UProbeAPI(
            protocol_config=Path(protocol),
            genomes_config=Path(genomes),
            output_dir=Path(output)
        )
        df_probes = pd.read_csv(probes)       
        df_processed = uprobe.post_process_probes(df_probes, raw_csv=raw, threads=threads)
        if df_processed.empty:
            log.warning("No probes remaining after post-processing!")
        else:
            log.info(f"Post-processing completed! {len(df_processed)} probes passed filters.")           
    except Exception as e:
        log.error(f"Post-processing failed: {e}")
        sys.exit(1)


def _generate_barcodes_dispatch(
    *,
    strategy: str,
    num_barcodes: int | None,
    length: int | None,
    alphabet: str,
    gc_limits: tuple[int, int] | None,
    prevent_patterns: list[str] | None,
    k_constraint: int | None,
) -> list:
    """Run barcode synthesis for the given CLI strategy (uses ``quick_generate`` where applicable)."""
    from .gen.barcodes import BarcodeGenerator, quick_generate

    ptn = prevent_patterns or None

    if strategy == "max_orthogonality":
        if not num_barcodes or not length:
            raise click.UsageError("'max_orthogonality' requires '--num-barcodes' and '--length'.")
        kw = {"alphabet": alphabet or "ACT", "rc_free": True}
        if gc_limits is not None:
            kw["gc_limits"] = gc_limits
        if ptn:
            kw["prevent_patterns"] = ptn
        return quick_generate(num_barcodes, length, **kw)

    if strategy == "pcr":
        if not num_barcodes:
            raise click.UsageError("'pcr' requires '--num-barcodes'.")
        length = length or 8
        lc, hc = length // 4, 3 * length // 4
        gc_pct = (round(100 * lc / length), round(100 * hc / length)) if length else (25, 75)
        return quick_generate(
            num_barcodes,
            length,
            alphabet="ACGT",
            rc_free=True,
            gc_limits=gc_pct,
            prevent_patterns=["AAAA", "TTTT", "CCCC", "GGGG"],
        )

    if strategy == "sequencing":
        if not num_barcodes:
            raise click.UsageError("'sequencing' requires '--num-barcodes'.")
        length = length or 12
        lc, hc = length // 3, 2 * length // 3
        gc_pct = (round(100 * lc / length), round(100 * hc / length)) if length else (33, 67)
        return quick_generate(
            num_barcodes,
            length,
            alphabet="ACGT",
            rc_free=True,
            gc_limits=gc_pct,
            prevent_patterns=["AAA", "TTT", "CCC", "GGG"],
        )

    if strategy == "max_size":
        if not k_constraint or not length:
            raise click.UsageError("'max_size' requires '--k-constraint' and '--length'.")
        gen = BarcodeGenerator(strategy="max_size")
        return gen.generate_max_size(
            length,
            k_constraint,
            alphabet=alphabet or "ACT",
            rc_free=True,
            gc_limits=None,
            prevent_patterns=ptn,
        )

    if strategy == "precomputed":
        raise click.UsageError("Strategy 'precomputed' is not implemented.")

    raise click.UsageError(f"Unknown strategy '{strategy}'.")


def _save_generated_barcodes(output_dir: Path, name: str, barcodes: list) -> None:
    csv_path = output_dir / f"{name}.csv"
    pd.DataFrame({"sequence": barcodes}).to_csv(csv_path, index=False)
    txt_path = output_dir / f"{name}.txt"
    txt_path.write_text("\n".join(str(b) for b in barcodes) + ("\n" if barcodes else ""), encoding="utf-8")


@cli.command(name='generate-barcodes')
@click.option('--protocol', '-p', type=click.Path(exists=True),
              help='Path to protocol config file (YAML). Used if no strategy is specified.')
@click.option('--output', '-o', default='./barcodes', type=click.Path(),
              help='Output directory for barcode files. [default: ./barcodes]')
@click.option('--strategy', '-s', type=click.Choice(['max_orthogonality', 'max_size', 'precomputed', 'pcr', 'sequencing']),
              help='Generation strategy to use.')
@click.option('--name', default='barcodes', help='Name for the generated barcode set.')
@click.option('--num-barcodes', type=int, help='Number of barcodes to generate.')
@click.option('--length', type=int, help='Length of each barcode.')
@click.option('--k-constraint', type=int, help='K-mer constraint for max_size strategy.')
@click.option('--library-name', help='Name of precomputed library (e.g., "kishi2018").')
@click.option('--alphabet', default='ACT', help='Nucleotide alphabet to use.')
@click.option('--gc-limits', help='GC content limits (min,max), e.g., "25,75".')
@click.option('--prevent-patterns', help='Comma-separated patterns to prevent, e.g., "AAAA,TTTT".')
@click.option('--save/--no-save', default=True, help='Save barcodes to a CSV file.')
@click.option('--analyze/--no-analyze', default=False, help='Analyze and save quality metrics.')
def generate_barcodes(protocol, output, strategy, name, num_barcodes, length, k_constraint,
                      library_name, alphabet, gc_limits, prevent_patterns, save, analyze):
    """
    Generate DNA barcode sequences.

    ``max_orthogonality``, ``pcr``, and ``sequencing`` use ``uprobe.core.gen.barcodes.quick_generate``
    (seqwalk ``max_orthogonality``). ``max_size`` uses ``BarcodeGenerator.generate_max_size``.

    With ``--save`` (default), writes ``{name}.csv`` (column ``sequence``) and ``{name}.txt``.
    Alternatively, pass ``-p/--protocol`` pointing to YAML that contains a ``barcode_generation`` block.
    """
    try:
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)

        gc_tuple: tuple[int, int] | None = None
        if gc_limits:
            parts = [int(x.strip()) for x in gc_limits.split(",")]
            if len(parts) != 2:
                raise click.UsageError('--gc-limits must be two comma-separated integers (percent GC), e.g. "25,75".')
            gc_tuple = (parts[0], parts[1])

        ptn_list: list[str] | None = None
        if prevent_patterns:
            ptn_list = [x.strip() for x in prevent_patterns.split(",") if x.strip()]

        strat = strategy
        nbc = num_barcodes
        ln = length
        kc = k_constraint
        ab = alphabet or "ACT"

        if protocol:
            log.info("Generating barcodes from protocol file...")
            cfg = yaml.safe_load(Path(protocol).read_text(encoding="utf-8"))
            bg = cfg.get("barcode_generation")
            if not bg:
                raise click.UsageError(
                    "Protocol YAML has no `barcode_generation` section. Add for example:\n\n"
                    "  barcode_generation:\n"
                    "    strategy: max_orthogonality\n"
                    "    num_barcodes: <N>\n"
                    "    length: <L>\n\n"
                    "Or omit `-p/--protocol` and pass `--strategy max_orthogonality` instead."
                )
            strat = bg.get("strategy") or strat or "max_orthogonality"
            nbc = bg.get("num_barcodes", nbc)
            ln = bg.get("length", ln)
            kc = bg.get("k_constraint", kc)
            ab = bg.get("alphabet", ab)
            glo = bg.get("gc_limits")
            if glo is not None:
                if isinstance(glo, (list, tuple)) and len(glo) == 2:
                    gc_tuple = (int(glo[0]), int(glo[1]))
                else:
                    raise click.UsageError("barcode_generation.gc_limits must be a two-element list.")
            plist = bg.get("prevent_patterns")
            if plist is not None:
                ptn_list = [str(x).strip() for x in plist if str(x).strip()]
        elif not strat:
            raise click.UsageError("Either `--strategy` or `-p/--protocol` (with `barcode_generation:`) must be provided.")

        _allowed = frozenset({"max_orthogonality", "max_size", "precomputed", "pcr", "sequencing"})
        if strat not in _allowed:
            raise click.UsageError(f"Unknown strategy {strat!r}. Choose from {sorted(_allowed)}.")

        if analyze:
            log.warning("'--analyze' is not implemented; skipping.")

        barcodes = _generate_barcodes_dispatch(
            strategy=strat,
            num_barcodes=nbc,
            length=ln,
            alphabet=ab,
            gc_limits=gc_tuple,
            prevent_patterns=ptn_list,
            k_constraint=kc,
        )
        if not barcodes:
            log.warning("No barcodes were generated.")
        elif save:
            _save_generated_barcodes(output_dir, name, barcodes)
            log.info(
                "Generated %s barcodes → %s, %s",
                len(barcodes),
                output_dir / f"{name}.csv",
                output_dir / f"{name}.txt",
            )
        else:
            log.info("Generated %s barcodes (--no-save); not writing files.", len(barcodes))
    except click.UsageError:
        raise
    except Exception as e:
        log.error(f"Barcode generation failed: {e}", exc_info=log.getEffectiveLevel() == logging.DEBUG)
        sys.exit(1)


@cli.command(name='generate-report')
@click.option('--protocol', '-p', required=True, type=click.Path(exists=True),
              help='Path to probe design protocol configuration file (YAML).')
@click.option('--genomes', '-g', required=True, type=click.Path(exists=True),
              help='Path to genome configuration file (YAML).')
@click.option('--probes', required=True, type=click.Path(exists=True),
              help='Path to processed probe results CSV file.')
@click.option('--output', '-o', default='./results', type=click.Path(),
              help='Output directory for report files. [default: ./results]')
@click.option('--no-plots', is_flag=True,
              help='Skip plot generation and only create text reports.')
@click.option('--pdf', is_flag=True, default=True,
              help='Generate PDF version of reports (default: enabled).')
@click.option('--no-pdf', is_flag=True,
              help='Skip PDF generation and only create markdown reports.')
def generate_report(protocol, genomes, probes, output, no_plots, pdf, no_pdf):
    """
    Generate interpretation report and plots for probe results.
    
    Creates detailed explanations of probe data columns and visualization plots
    to help users understand and select optimal probes.
    """
    try:
        import pandas as pd
        
        log.info("Generating probe analysis report...")
        uprobe = UProbeAPI(
            protocol_config=Path(protocol),
            genomes_config=Path(genomes),
            output_dir=Path(output)
        )
        df_probes = pd.read_csv(probes)       
        if df_probes.empty:
            log.error("No probe data found in the input file!")
            sys.exit(1)
        if pdf and not no_pdf:
            log.warning("PDF report generation is not supported by the current API; generating HTML only.")
        results = uprobe.generate_report(df_probes, include_plots=not no_plots)
        n_html = len(results.get('html_reports', []))
        if n_html == 0:
            log.warning("No reports or plots were generated. Check your protocol configuration.")
        else:
            log.info(f"Report generation completed!")
            log.info(f"Generated {n_html} HTML report(s)")
    except Exception as e:
        log.error(f"Report generation failed: {e}")
        sys.exit(1)


_DIVISION_CHOICE = click.Choice(['vertebrates', 'plants', 'fungi', 'metazoa', 'protists'])


def _dump_yaml(data) -> str:
    return yaml.safe_dump(data, sort_keys=False, default_flow_style=False, allow_unicode=True)


@cli.group(name='genome')
def genome():
    """
    Download reference genomes and inspect genomes.yaml.

    Sources: Ensembl (vertebrates, plus yeast, worm and fly) and Ensembl Genomes
    plants, fungi, metazoa and protists. Ensembl Bacteria is not supported.
    """


@genome.command(name='fetch')
@click.option('--species', '-s', required=True,
              help='Ensembl species name (homo_sapiens), display name (human) or taxonomy id (9606).')
@click.option('--release', '-r', default='latest',
              help='Ensembl release number, or "latest" (resolved and recorded). [default: latest]')
@click.option('--division', '-d', type=_DIVISION_CHOICE, default='vertebrates',
              help='Ensembl division; Ensembl Genomes divisions have their own release numbers. [default: vertebrates]')
@click.option('--out', '-o', default='./genomes', type=click.Path(file_okay=False),
              help='Directory for genome files; each genome gets a sub-directory. [default: ./genomes]')
@click.option('--genomes', '-g', default='genomes.yaml', type=click.Path(dir_okay=False),
              help='genomes.yaml to create or update. [default: genomes.yaml]')
@click.option('--name', help='Entry name in genomes.yaml. [default: <species>_<assembly>_e<release>]')
@click.option('--align-index', help='Comma-separated aligners for the entry, e.g. "bowtie2,blast". '
                                    '[default: keep the existing value, else bowtie2, plus bowtie if '
                                    'bowtie-build (Bowtie 1) is installed]')
@click.option('--jellyfish/--no-jellyfish', default=None,
              help='Enable kmer_count for the entry. [default: keep the existing value, else on if jellyfish is installed]')
@click.option('--build-index', is_flag=True, help='Build aligner indices afterwards (same as `uprobe build-index`).')
@click.option('--threads', '-t', default=10, type=int, help='Threads for --build-index. [default: 10]')
@click.option('--keep-compressed', is_flag=True, help='Keep the downloaded .gz files.')
@click.option('--force', is_flag=True,
              help='Re-download verified files, and replace a different genome already in the target directory.')
@click.option('--dry-run', is_flag=True, help='Resolve release, species and URLs, print the plan, download nothing.')
@click.option('--mirror', help='Base URL replacing https://ftp.ensembl.org/pub '
                               '(or https://ftp.ensemblgenomes.ebi.ac.uk/pub for other divisions).')
def genome_fetch(species, release, division, out, genomes, name, align_index, jellyfish, build_index,
                 threads, keep_compressed, force, dry_run, mirror):
    """
    Download an Ensembl genome FASTA and GTF and register them in genomes.yaml.

    Prefers the primary assembly FASTA and falls back to toplevel. Files are
    verified against Ensembl's CHECKSUMS, decompressed, and recorded with their
    release and URLs. Re-running skips verified files and resumes partial
    downloads. The resulting entry is printed to stdout as YAML.

    Example: uprobe genome fetch --species homo_sapiens --release latest --build-index
    """
    from .genome_fetch import GenomeFetchError, fetch_genome
    aligners = [a.strip() for a in align_index.split(',') if a.strip()] if align_index else None
    try:
        result = fetch_genome(
            species, release=release, division=division, out_dir=out, genomes_yaml=genomes, name=name,
            align_index=aligners, jellyfish=jellyfish, build_index=build_index, threads=threads,
            keep_compressed=keep_compressed, force=force, dry_run=dry_run, mirror=mirror,
        )
    except GenomeFetchError as e:
        log.error(str(e))
        sys.exit(1)
    except Exception as e:
        log.error(f"Genome fetch failed: {e}")
        sys.exit(1)
    if dry_run:
        click.echo("# Dry run: nothing was downloaded.")
        click.echo(_dump_yaml({'plan': result['plan'].to_dict(), 'directory': str(result['directory']),
                               'genomes_yaml': str(result['genomes_yaml']),
                               'entry': {result['name']: result['entry']}}), nl=False)
        return
    click.echo(f"# Registered in {result['genomes_yaml']}; use `genome: {result['name']}` in your protocol.")
    if not build_index:
        click.echo(f"# Build indices with: uprobe build-index -p <protocol.yaml> -g {genomes} "
                   "(or re-run with --build-index).")
    click.echo(_dump_yaml({result['name']: result['entry']}), nl=False)


@genome.command(name='list')
@click.option('--genomes', '-g', default='genomes.yaml', type=click.Path(exists=True, dir_okay=False),
              help='Path to genome configuration file (YAML). [default: genomes.yaml]')
@click.option('--json', 'as_json', is_flag=True, help='Print machine-readable JSON.')
def genome_list(genomes, as_json):
    """Show genomes.yaml entries and whether their files and indices exist."""
    import json
    from .genome_fetch import GenomeFetchError, describe_genomes
    try:
        rows = describe_genomes(genomes)
    except (GenomeFetchError, yaml.YAMLError) as e:
        log.error(f"Cannot read {genomes}: {e}")
        sys.exit(1)
    if as_json:
        click.echo(json.dumps(rows, indent=2))
        return
    if not rows:
        click.echo(f"No genomes in {genomes}.")
    state = {True: 'ready', False: 'missing', None: 'unknown'}
    for row in rows:
        click.echo(row['name'])
        if row['description']:
            click.echo(f"  description: {row['description']}")
        click.echo(f"  fasta: {row['fasta']} ({'ok' if row['fasta_exists'] else 'MISSING'})")
        click.echo(f"  gtf:   {row['gtf']} ({'ok' if row['gtf_exists'] else 'MISSING'})")
        for aligner, status in row['indexes'].items():
            click.echo(f"  {aligner} index: " + ", ".join(f"{kind} {state[value]}" for kind, value in status.items()))
        click.echo(f"  jellyfish: {'on' if row['jellyfish'] else 'off'}")
        source = row['source'] or {}
        if source.get('release'):
            click.echo(f"  source: {source.get('provider', 'Ensembl')} {source.get('division', '')} "
                       f"release {source['release']} ({source.get('fasta_url')})")


@genome.command(name='search')
@click.argument('query')
@click.option('--division', '-d', type=_DIVISION_CHOICE, default='vertebrates',
              help='Ensembl division to search. [default: vertebrates]')
@click.option('--release', '-r', default='latest', help='Ensembl release number or "latest". [default: latest]')
@click.option('--limit', default=20, type=int, help='Maximum number of matches. [default: 20]')
@click.option('--mirror', help='Base URL replacing the Ensembl FTP "pub" URL.')
def genome_search(query, division, release, limit, mirror):
    """Find Ensembl species names by name, common name or taxonomy id."""
    from .genome_fetch import GenomeFetchError, search
    try:
        release, rows = search(query, division=division, release=release, mirror=mirror, limit=limit)
    except GenomeFetchError as e:
        log.error(str(e))
        sys.exit(1)
    if not rows:
        click.echo(f"No species matching '{query}' in Ensembl {division} release {release}.")
        sys.exit(1)
    for row in rows:
        click.echo(f"{row.get('species')}\t{row.get('name', '')}\t{row.get('assembly', '')}\t"
                   f"taxid {row.get('taxonomy_id', '')}")


@cli.group()
def templates():
    """
    Protocol templates for published probe technologies.

    Each template is a complete protocol, including the source protocol's
    hybridization conditions; its `meta` block lists the cited sources and
    what to replace.
    """


@templates.command(name='list')
def templates_list():
    """List available templates with technology and primary citation."""
    from uprobe.templates import list_templates
    rows = list_templates()
    headers = {'name': 'NAME', 'technology': 'TECHNOLOGY', 'citation': 'CITATION'}
    widths = {k: max(len(r[k]) for r in [headers, *rows]) for k in headers}
    for row in [headers, *rows]:
        click.echo(f"{row['name']:<{widths['name']}}  {row['technology']:<{widths['technology']}}  {row['citation']}")


@templates.command(name='show')
@click.argument('name')
@click.option('--output', '-o', type=click.Path(dir_okay=False),
              help='Write the template to this file instead of stdout.')
@click.option('--force', is_flag=True, help='Overwrite the output file if it exists.')
def templates_show(name, output, force):
    """Print template NAME, or save it with -o as a starting protocol."""
    from uprobe.templates import read_template
    try:
        text = read_template(name)
    except KeyError as e:
        raise click.BadParameter(e.args[0], param_hint='NAME')
    if not output:
        click.echo(text, nl=False)
        return
    path = Path(output)
    if path.exists() and not force:
        raise click.UsageError(f"{path} exists; use --force to overwrite.")
    path.write_text(text, encoding='utf-8')
    click.echo(f"Template {name!r} written to {path}", err=True)


@cli.group()
def skill():
    """
    The U-Probe Agent Skill for AI coding agents (Codex, Claude Code, ...).

    The skill (SKILL.md, references/, scripts/) teaches an agent the U-Probe
    workflow and the probe design rules behind it. `uprobe skill install`
    copies it to where an agent discovers skills.
    """


@skill.command(name='path')
def skill_path_command():
    """Print the skill folder shipped with this U-Probe installation."""
    from uprobe.skills import skill_path
    click.echo(str(skill_path()))


@skill.command(name='install')
@click.option('--target', 'targets', multiple=True, metavar='codex|claude|DIR',
              help='Where to install: "codex" (.agents/skills), "claude" (.claude/skills), or a skills '
                   'directory DIR (the skill goes to DIR/u-probe). Repeatable. [default: codex and claude]')
@click.option('--project/--user', 'project', default=True,
              help='For codex and claude: the project directory (the working directory) or the user\'s '
                   'home directory (~/.agents/skills, ~/.claude/skills). [default: --project]')
@click.option('--force', is_flag=True, help='Replace an existing u-probe skill folder (e.g. after upgrading U-Probe).')
def skill_install(targets, project, force):
    """
    Copy the U-Probe skill to where coding agents discover skills.

    \b
    Codex:       .agents/skills/u-probe   (project)   ~/.agents/skills/u-probe   (--user)
    Claude Code: .claude/skills/u-probe   (project)   ~/.claude/skills/u-probe   (--user)

    Project installs go under the working directory; run it at the repository root.
    Run again with --force after upgrading U-Probe.

    Example: uprobe skill install --target claude --user
    """
    from uprobe.skills import SkillInstallError, install_skill, skills_dir
    targets = targets or ('codex', 'claude')
    try:
        directories = [skills_dir(target, user=not project) for target in targets]
    except SkillInstallError as e:
        raise click.UsageError(str(e))
    for target, directory in zip(targets, directories):
        try:
            dest = install_skill(directory, force=force)
        except (SkillInstallError, OSError) as e:
            raise click.ClickException(str(e))
        label = target if target in ('codex', 'claude') else 'directory'
        click.echo(f"Installed the u-probe skill (U-Probe {__version__}) for {label}: {dest}")


@cli.command()
def version():
    """Show the U-Probe version."""
    click.echo(f"U-Probe version {__version__}")


@cli.command()
@click.option('--host', default=None, help='Host to bind the server to (overrides .env).')
@click.option('--port', default=None, type=int, help='Port to bind the server to (overrides .env).')
@click.option('--workers', default=None, type=int, help='Number of worker processes (overrides .env).')
@click.option('--env', default=None, type=click.Choice(['development', 'production']), help='Environment mode (overrides APP_ENV).')
def server(host, port, workers, env):
    """Start the U-Probe HTTP web server using Uvicorn."""
    try:
        import os
        import uvicorn
        
        # Override environment variables if CLI arguments are provided
        if env is not None:
            os.environ["APP_ENV"] = env
        if host is not None:
            os.environ["HOST"] = host
        if port is not None:
            os.environ["PORT"] = str(port)
        if workers is not None:
            os.environ["WORKERS"] = str(workers)
            
        # Import start_server from server module to reuse the logic
        from uprobe.http.server import start_server
        
        start_server()
        
    except ImportError as e:
        log.error(f"Failed to start server: {e}")
        log.error("Please ensure web dependencies (fastapi, uvicorn) are installed.")
        sys.exit(1)


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()
