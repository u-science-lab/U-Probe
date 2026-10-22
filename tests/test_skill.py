"""The Agent Skill in uprobe/skills/u-probe: SKILL.md and its routing to the references, the
installer (`uprobe skill path`, `uprobe skill install`) and scripts/check_probes.py."""
import importlib.util
import json
import os
import random
import re
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from uprobe.core.cli import cli
from uprobe.core.gen.fun import extract_trans_seqs
from uprobe.core.utils import reverse_complement
from uprobe.skills import SKILL_NAME, install_skill, read_frontmatter, skill_path, skills_dir

SKILL = skill_path()
SKILL_MD = SKILL / 'SKILL.md'
# Agent Skills specification (agentskills.io/specification): name and description rules.
NAME_RULE = re.compile(r'^[a-z0-9]+(-[a-z0-9]+)*$')
# Files of the skill named in a text: references/... .md and scripts/... .py.
NAMED_FILE = re.compile(r'(?:references|scripts)/[\w./-]+?\.(?:md|py)')

spec = importlib.util.spec_from_file_location('check_probes', SKILL / 'scripts' / 'check_probes.py')
check_probes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_probes)


def files_of(root):
    return sorted(p.relative_to(root).as_posix() for p in root.rglob('*')
                  if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc')


# -- SKILL.md -----------------------------------------------------------------------------

def test_frontmatter_follows_the_agent_skills_spec():
    assert SKILL_MD.read_text(encoding='utf-8').startswith('---\n')
    meta = read_frontmatter(SKILL_MD)
    assert meta['name'] == SKILL_NAME == SKILL.name
    assert NAME_RULE.match(meta['name']) and len(meta['name']) <= 64
    description = meta['description']
    assert isinstance(description, str) and 0 < len(description) <= 1024
    for trigger in ('smFISH', 'HCR', 'MERFISH', 'seqFISH', 'padlock', 'ISS', 'STARmap', 'Oligopaint',
                    'DNA FISH', 'oligo pools'):
        assert trigger in description, trigger


def test_skill_md_stays_short():
    # The body is loaded whenever the skill triggers; the references carry the detail.
    assert len(SKILL_MD.read_text(encoding='utf-8').splitlines()) <= 260


def test_every_file_named_in_skill_md_exists():
    named = set(NAMED_FILE.findall(SKILL_MD.read_text(encoding='utf-8')))
    assert len(named) >= 19
    assert sorted(p for p in named if not (SKILL / p).is_file()) == []


def test_every_reference_is_routed_from_skill_md():
    text = SKILL_MD.read_text(encoding='utf-8')
    index = (SKILL / 'references' / 'techniques' / 'README.md').read_text(encoding='utf-8')
    for path in sorted((SKILL / 'references').rglob('*.md')):
        relative = path.relative_to(SKILL).as_posix()
        # A technique file may be reached through the techniques index.
        technique = path.parent.name == 'techniques' and path.name != 'README.md'
        assert relative in text or (technique and f'`{path.name}`' in index), relative


def test_skill_md_names_only_existing_commands_and_templates():
    from uprobe.templates import template_names
    text = SKILL_MD.read_text(encoding='utf-8')
    for name in re.findall(r'uprobe templates show (\w+)', text):
        assert name in template_names()
    commands = set(cli.commands)
    for command in re.findall(r'(?:^|`)uprobe ([a-z][a-z-]*)', text, flags=re.MULTILINE):
        assert command in commands, command


# -- What the skill folder holds -------------------------------------------------------------
# The folder ships in the wheel and `uprobe skill install` copies all of it, so an agent sees every
# file in it. Notes, corpora and instructions for other agents belong outside the package (e.g.
# fish-methods-corpus/ at the repository root).

# Left out by the installer (uprobe.skills._IGNORE) and by the wheel.
NOT_SHIPPED = {'__pycache__', '.DS_Store'}
AGENT_INSTRUCTIONS = {'agents.md', 'claude.md'}
BINARY_SUFFIXES = {'.pdf', '.xlsx', '.xls', '.docx', '.doc', '.pptx', '.ppt', '.zip', '.gz', '.png', '.jpg',
                   '.jpeg', '.gif', '.h5', '.npz', '.pkl'}


def shipped_files():
    return sorted(p for p in SKILL.rglob('*') if p.is_file() and not NOT_SHIPPED & set(p.parts)
                  and p.suffix != '.pyc')


def test_skill_folder_holds_only_skill_md_references_and_scripts():
    top = sorted(p.name for p in SKILL.iterdir() if p.name not in NOT_SHIPPED)
    assert top == ['SKILL.md', 'references', 'scripts']
    assert SKILL_MD.is_file() and (SKILL / 'references').is_dir() and (SKILL / 'scripts').is_dir()
    # MANIFEST.in ships only references/**/*.md and scripts/**/*.py; anything else would be
    # installed from a checkout but missing from the wheel.
    stray = [p.relative_to(SKILL).as_posix() for p in shipped_files()
             if p != SKILL_MD and not (p.parts[len(SKILL.parts)] == 'references' and p.suffix == '.md')
             and not (p.parts[len(SKILL.parts)] == 'scripts' and p.suffix == '.py')]
    assert stray == []


def test_skill_folder_has_no_agent_instructions_or_binaries():
    found = []
    for path in shipped_files():
        binary = path.suffix.lower() in BINARY_SUFFIXES or b'\0' in path.read_bytes()[:8192]
        if path.name.lower() in AGENT_INSTRUCTIONS or binary:
            found.append(path.relative_to(SKILL).as_posix())
    assert found == []


# -- uprobe skill path / install ------------------------------------------------------------

def test_skill_path_prints_the_packaged_folder():
    result = CliRunner().invoke(cli, ['skill', 'path'])
    assert result.exit_code == 0, result.output
    assert Path(result.output.strip()) == SKILL and SKILL_MD.is_file()


def test_install_into_a_directory_gives_a_valid_skill(tmp_path):
    skills = tmp_path / 'skills'
    result = CliRunner().invoke(cli, ['skill', 'install', '--target', str(skills)])
    assert result.exit_code == 0, result.output
    dest = skills / SKILL_NAME
    assert str(dest) in result.output
    assert read_frontmatter(dest / 'SKILL.md') == read_frontmatter(SKILL_MD)
    assert files_of(dest) == files_of(SKILL)
    named = set(NAMED_FILE.findall((dest / 'SKILL.md').read_text(encoding='utf-8')))
    assert named and all((dest / p).is_file() for p in named)
    assert list(skills.iterdir()) == [dest]  # no staging folder left behind


def test_install_replaces_only_with_force(tmp_path):
    runner = CliRunner()
    assert runner.invoke(cli, ['skill', 'install', '--target', str(tmp_path)]).exit_code == 0
    (tmp_path / SKILL_NAME / 'stale.md').write_text('from an older version')
    again = runner.invoke(cli, ['skill', 'install', '--target', str(tmp_path)])
    assert again.exit_code == 1 and 'use --force' in again.output
    assert (tmp_path / SKILL_NAME / 'stale.md').exists()
    forced = runner.invoke(cli, ['skill', 'install', '--target', str(tmp_path), '--force'])
    assert forced.exit_code == 0, forced.output
    assert files_of(tmp_path / SKILL_NAME) == files_of(SKILL)


def test_install_never_replaces_a_folder_that_is_not_the_skill(tmp_path):
    other = tmp_path / SKILL_NAME
    other.mkdir()
    (other / 'notes.txt').write_text('mine')
    result = CliRunner().invoke(cli, ['skill', 'install', '--target', str(tmp_path), '--force'])
    assert result.exit_code == 1 and 'not a U-Probe skill' in result.output
    assert files_of(other) == ['notes.txt']


@pytest.mark.skipif(os.name == 'nt', reason='symbolic links need privileges on Windows')
def test_install_replaces_a_symlink_without_following_it(tmp_path):
    checkout = tmp_path / 'checkout'
    checkout.mkdir()
    (checkout / 'keep.txt').write_text('x')
    skills = tmp_path / 'skills'
    skills.mkdir()
    (skills / SKILL_NAME).symlink_to(checkout, target_is_directory=True)
    dest = install_skill(skills, force=True)
    assert not dest.is_symlink() and (dest / 'SKILL.md').is_file()
    assert files_of(checkout) == ['keep.txt']


def test_skills_dirs_of_the_agents(tmp_path):
    # Codex: developers.openai.com/codex/skills; Claude Code: code.claude.com/docs/en/skills.
    project, home = tmp_path / 'project', tmp_path / 'home'
    assert skills_dir('codex', project_dir=project) == project / '.agents' / 'skills'
    assert skills_dir('claude', project_dir=project) == project / '.claude' / 'skills'
    assert skills_dir('codex', user=True, home=home) == home / '.agents' / 'skills'
    assert skills_dir('claude', user=True, home=home) == home / '.claude' / 'skills'
    assert skills_dir(str(tmp_path / 'x')) == tmp_path / 'x'


def test_install_for_codex_and_claude(tmp_path, monkeypatch):
    project, home = tmp_path / 'project', tmp_path / 'home'
    project.mkdir()
    home.mkdir()
    monkeypatch.chdir(project)
    monkeypatch.setenv('HOME', str(home))
    monkeypatch.setenv('USERPROFILE', str(home))
    runner = CliRunner()
    both = runner.invoke(cli, ['skill', 'install'])  # default: both agents, project level
    assert both.exit_code == 0, both.output
    assert (project / '.agents/skills/u-probe/SKILL.md').is_file()
    assert (project / '.claude/skills/u-probe/SKILL.md').is_file()
    for target, folder in (('codex', '.agents'), ('claude', '.claude')):
        result = runner.invoke(cli, ['skill', 'install', '--target', target, '--user'])
        assert result.exit_code == 0, result.output
        assert (home / folder / 'skills/u-probe/SKILL.md').is_file()
    bad = runner.invoke(cli, ['skill', 'install', '--target', str(tmp_path / 'x'), '--user'])
    assert bad.exit_code == 2 and '--user applies' in bad.output


def _check_commands(text):
    """The check_probes.py command lines of the fenced bash blocks of a reference, continuations joined."""
    commands = []
    for block in re.findall(r'```bash\n(.*?)```', text, flags=re.S):
        joined = re.sub(r'\\\n\s*', ' ', block)
        commands += [line.strip() for line in joined.splitlines() if 'scripts/check_probes.py' in line]
    return commands


TECHNIQUES = sorted(p for p in (SKILL / 'references' / 'techniques').glob('*.md') if p.name != 'README.md')


@pytest.mark.parametrize('path', TECHNIQUES, ids=lambda p: p.name)
def test_every_technique_checklist_ends_with_a_check_probes_command(path, tmp_path):
    text = path.read_text(encoding='utf-8')
    checklist = re.split(r'^## \d+\. (?:Verification checklist|Checklist)\s*$', text, flags=re.M)
    assert len(checklist) == 2, f'{path.name}: no checklist section'
    commands = _check_commands(checklist[1].split('\n## ')[0])
    assert commands, f'{path.name}: the checklist has no check_probes.py command'
    for command in commands:  # the options must parse with the script's own argument parser
        tokens = command.replace('<genome dir>', 'genomes').split()
        assert tokens[:2] == ['python', '<skill>/scripts/check_probes.py'], command
        args = check_probes.parse_args(tokens[2:])
        assert args.tsv.endswith('.tsv'), command


def test_the_standard_length_recipe_of_troubleshooting_validates_without_warnings(tmp_path):
    from uprobe.core.validate import validate_files
    text = (SKILL / 'references' / 'uprobe' / 'troubleshooting.md').read_text(encoding='utf-8')
    blocks = [b for b in re.findall(r'```yaml\n(.*?)```', text, flags=re.S) if 'prefer: len_off' in b]
    assert len(blocks) == 1
    path = tmp_path / 'standard_length_first.yaml'
    path.write_text(blocks[0], encoding='utf-8')
    assert [f.message for f in validate_files(path) if f.severity in ('error', 'warning')] == []


# -- scripts/check_probes.py on hand-made tables ---------------------------------------------

_rng = random.Random(3)
RNA = ''.join(_rng.choice('ACGT') for _ in range(400))
COLUMNS = ['probe_id', 'target', 'target_start', 'target_end', 'target_region', 'probe', 'probe.binding',
           'probe_tm', 'target_mapped_genes', 'target_mapped_genes_off_targets']


def rows_at(starts, length=20, target='G'):
    rows = []
    for n, start in enumerate(starts):
        window = RNA[start:start + length]
        rows.append({'probe_id': f'{target}_{n}', 'target': target, 'target_start': start + 1,
                     'target_end': start + length, 'target_region': window, 'probe': reverse_complement(window),
                     'probe.binding': reverse_complement(window), 'probe_tm': 55.0, 'target_mapped_genes': 1,
                     'target_mapped_genes_off_targets': ''})
    return rows


def run_check(tmp_path, rows, *options, columns=COLUMNS):
    table = tmp_path / 'probes.tsv'
    table.write_text('\n'.join(['\t'.join(columns)] + ['\t'.join(str(r.get(c, '')) for c in columns)
                                                         for r in rows]) + '\n')
    return check_probes.check(check_probes.parse_args([str(table), *options]))


def test_check_passes_reverse_complement_probes(tmp_path):
    result = run_check(tmp_path, rows_at([0, 22, 44, 70]))
    assert result['ok'] and result['failures'] == [] and result['warnings'] == []
    entry = result['targets']['G']
    assert entry['probes'] == 4 and entry['gap'] == {'min': 2, 'median': 2, 'max': 6} and entry['overlaps'] == 0
    assert entry['parts']['probe.binding']['length'] == {'min': 20, 'median': 20, 'max': 20}
    assert entry['tm']['probe_tm']['median'] == 55.0


def test_check_fails_sense_probes_and_placeholders(tmp_path):
    rows = rows_at([0, 30, 60])
    rows[0]['probe.binding'] = rows[0]['probe'] = rows[0]['target_region']  # the mRNA itself
    rows[1]['probe'] = rows[1]['probe'][:-3] + 'NNN'
    result = run_check(tmp_path, rows)
    assert not result['ok']
    assert any('G_0' in f and 'it is the sense sequence' in f for f in result['failures'])
    assert any(f.startswith('G_1: probe is not pure ACGT') for f in result['failures'])
    # Where the letter is: the first 30 nt of a long probe need not contain it.
    assert "G_1: probe is not pure ACGT: 'N' at position 18 of 20" in result['failures']


def test_check_limits_and_off_targets(tmp_path):
    rows = rows_at([0, 15, 40])  # windows 1-20 and 16-35 overlap
    rows[2].update(target_mapped_genes=2, target_mapped_genes_off_targets='OTHER:1')
    result = run_check(tmp_path, rows)
    assert result['ok']
    assert any('1 overlapping windows' in w for w in result['warnings'])
    assert any('OTHER:1' in w for w in result['warnings'])
    assert result['targets']['G']['off_targets']['target_mapped_genes'] == {'max_genes': 2, 'rows_with_off_targets': 1}
    strict = run_check(tmp_path, rows, '--min-gap', '2', '--tm', '56', '60', '--gc', '0.9', '1', '--min-probes', '5')
    assert not strict['ok'] and len(strict['failures']) == 4


def test_check_sense_parts_of_cdna_padlocks(tmp_path):
    rows = rows_at([0, 40])
    for row in rows:  # HybISS: arm_5 = window[15:30], arm_3 = window[0:15], both mRNA sense
        row['padlock.arm_5'], row['padlock.arm_3'] = row['target_region'][15:], row['target_region'][:15]
    columns = COLUMNS[:5] + ['padlock.arm_5', 'padlock.arm_3']
    ok = run_check(tmp_path, rows, '--sense', 'padlock.arm_5', '--sense', 'padlock.arm_3', columns=columns)
    assert ok['ok'], ok['failures']
    wrong = run_check(tmp_path, rows, '--binding', 'padlock.arm_5', columns=columns)
    assert not wrong['ok']


def test_check_against_a_transcript_fasta(tmp_path):
    fasta = tmp_path / 'tx.fa'
    fasta.write_text(f'>ENSG1_ENST1\n{RNA[:60]}\n{RNA[60:300]}\n')
    rows = rows_at([0, 100, 350])  # the last window lies beyond the transcript
    result = run_check(tmp_path, rows, '--transcripts', str(fasta))
    assert [f for f in result['failures'] if 'not on the designed transcript' in f] == \
        ['G_2: window not on the designed transcript (ENSG1_ENST1)']
    assert result['targets']['G']['span'] == {'ENSG1_ENST1': {'from': 1, 'to': 120, 'length': 300, 'fraction': 0.4}}
    named = run_check(tmp_path, rows[:2], '--transcripts', str(fasta), '--transcript', 'G=ENST1.4')
    assert named['ok'] and set(named['targets']['G']['span']) == {'ENST1'}


def test_check_spacing_of_genomic_windows(tmp_path):
    # source: genome tables (Oligopaint) place windows by sub_region "s-e" along the region and
    # have no target_start: --min-gap passed without checking anything.
    columns = ['probe_id', 'target', 'sub_region', 'target_region', 'oligo.homology']
    rows = [{'probe_id': f'chr1:1-400_{n}', 'target': 'chr1:1-400', 'sub_region': f'{s + 1}-{s + 35}',
             'target_region': RNA[s:s + 35], 'oligo.homology': RNA[s:s + 35]} for n, s in enumerate([0, 30, 100])]
    result = run_check(tmp_path, rows, '--sense', 'oligo.homology', '--min-gap', '0', columns=columns)
    assert result['targets']['chr1:1-400']['gap'] == {'min': -5, 'median': 15, 'max': 35}
    assert result['failures'] == ['chr1:1-400: smallest gap -5 nt < --min-gap 0']
    unplaced = run_check(tmp_path, rows, '--sense', 'oligo.homology', '--min-gap', '0', columns=columns[:2] + columns[3:])
    assert not unplaced['ok'] and 'cannot check --min-gap' in unplaced['failures'][0]


def test_check_per_part_gc_and_tm(tmp_path):
    # A split unit: a window-level GC passes while one half is AT-rich (SNAIL-style padlock + primer).
    at, gc_rich = 'ATATTAATTATAATTATATA', 'GCGGCCGCGCGGCCGCGCGC'
    rows = []
    for n in range(2):
        rows.append({'probe_id': f'U_{n}', 'target': 'G', 'target_start': 1 + 40 * n, 'target_end': 40 + 40 * n,
                     'target_region': at + gc_rich, 'half_a': reverse_complement(at), 'half_b': reverse_complement(gc_rich),
                     'half_a_tm': 41.5, 'half_b_tm': 66.0})
    columns = ['probe_id', 'target', 'target_start', 'target_end', 'target_region', 'half_a', 'half_b',
               'half_a_tm', 'half_b_tm']
    common = ['--binding', 'half_a', '--binding', 'half_b']
    whole = run_check(tmp_path, rows, *common, '--gc', '0.3', '1', columns=columns)
    assert not whole['ok'] and any('half_a outside GC' in f for f in whole['failures'])  # --gc applies to every part
    ok = run_check(tmp_path, rows, *common, '--part-gc', 'half_b', '0.3', '1', '--part-tm', 'half_b_tm', '55', '70',
                   columns=columns)
    assert ok['ok'] and ok['targets']['G']['part_gc']['half_b']['median'] == 1.0
    bad = run_check(tmp_path, rows, *common, '--part-gc', 'half_a', '0.3', '0.6', '--part-tm', 'half_a_tm', '50', '70',
                    columns=columns)
    assert sorted(bad['failures']) == ['G: 2 half_a outside GC 0.3-0.6', 'G: 2 half_a_tm outside Tm 50-70']
    with pytest.raises(SystemExit, match='no column'):
        run_check(tmp_path, rows, *common, '--part-gc', 'nope', '0.3', '0.6', columns=columns)


def test_check_reports_duplicates_and_multi_site_windows(tmp_path):
    # Tandem copies: the same 20-mer occurs three times on the transcript but is one oligo.
    unit = RNA[:20]
    transcript = RNA[100:130] + unit + RNA[130:150] + unit + RNA[150:170] + unit + RNA[170:200]
    fasta = tmp_path / 'tx.fa'
    fasta.write_text(f'>ENSG1_ENST1\n{transcript}\n')
    other = transcript[0:20]
    rows = []
    for n, window in enumerate([unit, unit, other]):  # the table lists the repeated oligo twice
        start = transcript.find(window)
        rows.append({'probe_id': f'G_{n}', 'target': 'G', 'target_start': start + 1 + 100 * (n == 1),
                     'target_end': start + 20 + 100 * (n == 1), 'target_region': window,
                     'probe': reverse_complement(window), 'probe.binding': reverse_complement(window)})
    columns = ['probe_id', 'target', 'target_start', 'target_end', 'target_region', 'probe', 'probe.binding']
    result = run_check(tmp_path, rows, '--transcripts', str(fasta), '--transcript', 'G=ENST1', '--min-gap', '12',
                       '--min-probes', '3', columns=columns)
    entry = result['targets']['G']
    assert entry['probes'] == 3 and entry['distinct_oligos'] == 2 and entry['duplicate_rows'] == 1
    assert entry['multi_site'][0]['sites'] == 3 and len(entry['multi_site']) == 1
    assert any('2 distinct oligos' in w for w in result['warnings'])
    assert any('more than one place' in w for w in result['warnings'])
    assert any('fewer than --min-probes 3' in f and '2 distinct' in f for f in result['failures'])
    # --min-gap also looks at every site: the other window lies 9 nt before the first copy (gap 9 < 12).
    assert any('over all their sites' in f for f in result['failures'])
    # Without duplicates or repeats nothing is reported.
    clean = run_check(tmp_path, rows_at([0, 30]))
    assert 'multi_site' not in clean['targets']['G'] and clean['targets']['G']['distinct_oligos'] == 2


def test_check_command_line(tmp_path, capsys):
    table = tmp_path / 'probes.tsv'
    rows = rows_at([0, 30])
    table.write_text('\n'.join(['\t'.join(COLUMNS)] + ['\t'.join(str(r[c]) for c in COLUMNS) for r in rows]) + '\n')
    assert check_probes.main([str(table), '--json']) == 0
    assert json.loads(capsys.readouterr().out)['targets']['G']['probes'] == 2
    assert check_probes.main([str(table)]) == 0
    assert 'OK (0 failures, 0 warnings)' in capsys.readouterr().out
    (tmp_path / 'other.tsv').write_text('a\tb\n1\t2\n')
    assert check_probes.main([str(tmp_path / 'other.tsv')]) == 2


# -- scripts/check_probes.py on a real U-Probe run -------------------------------------------

MANE = 'tag "MANE_Select"; tag "Ensembl_canonical";'
CANONICAL = 'tag "Ensembl_canonical";'
# (gene, transcript, tags, start, end, strand): P's ENSTP2 has a long alternative last exon,
# sequence its MANE Select transcript ENSTP1 lacks; M is on the - strand.
EXONS = [('P', 'ENSTP1', MANE, 101, 400, '+'), ('P', 'ENSTP1', MANE, 501, 800, '+'),
         ('P', 'ENSTP2', '', 101, 400, '+'), ('P', 'ENSTP2', '', 1001, 2000, '+'),
         ('M', 'ENSTM3', CANONICAL, 2101, 2300, '-'), ('M', 'ENSTM3', CANONICAL, 2401, 2600, '-')]


@pytest.fixture
def toy(tmp_path):
    rng = random.Random(11)
    genome = ''.join(rng.choice('ACGT') for _ in range(2700))
    (tmp_path / 'g.fa').write_text(f'>1\n{genome}\n')
    (tmp_path / 'g.gtf').write_text(''.join(
        f'1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "G{gene}"; transcript_id "{transcript}"; '
        f'gene_name "{gene}"; {tags}\n' for gene, transcript, tags, start, end, strand in EXONS))
    (tmp_path / 'genomes.yaml').write_text(yaml.safe_dump(
        {'toy': {'fasta': str(tmp_path / 'g.fa'), 'gtf': str(tmp_path / 'g.gtf'), 'align_index': []}}))
    extract_trans_seqs(str(tmp_path / 'g.gtf'), str(tmp_path / 'g.fa'), str(tmp_path / 'transcript.fa'))
    return tmp_path


def design(toy, transcripts):
    protocol = {
        'name': 'smfish', 'genome': 'toy', 'targets': ['P', 'M'], 'encoding': {},
        'extracts': {'target_region': {'source': 'mRNA', 'transcripts': transcripts, 'length': 20, 'step': 1}},
        'probes': {'smfish_probe': {'template': '{binding}', 'parts': {'binding': {'expr': 'rc(target_region)'}}}},
        'hybridization': {'temperature_c': 37, 'na_mM': 390, 'formamide_pct': 10},
        'attributes': {'target_gc': {'target': 'target_region', 'type': 'gc_content'},
                       'binding_tm': {'target': 'smfish_probe.binding', 'type': 'annealing_temperature'}},
        'post_process': {'remove_overlap': {'location_interval': 2}, 'equal_space': {'number_desired': 12}},
        'summary': {'report_name': 'rna_report', 'attributes': ['target_gc', 'binding_tm']}}
    path = toy / f'{transcripts}.yaml'
    path.write_text(yaml.safe_dump(protocol))
    output = toy / f'out_{transcripts}'
    result = CliRunner().invoke(cli, ['run', '-p', str(path), '-g', str(toy / 'genomes.yaml'),
                                      '-o', str(output), '--json'])
    assert result.exit_code == 0, result.output
    summary = json.loads((output / 'summary.json').read_text())
    return output / summary['output_files']['tsv'], output / 'summary.json'


def test_check_a_run_on_the_reference_transcripts(toy, capsys):
    table, summary = design(toy, 'mane')
    code = check_probes.main([str(table), '--summary', str(summary), '--transcripts', str(toy / 'transcript.fa'),
                              '--min-gap', '2', '--json'])
    result = json.loads(capsys.readouterr().out)
    assert code == 0 and result['ok'], result['failures']
    assert result['binding_columns'] == ['smfish_probe.binding']
    assert result['targets']['P']['transcripts'] == {'ENSTP1': 600}
    assert result['targets']['P']['chosen_by'] == 'MANE_Select' and result['targets']['P']['probes'] == 12
    assert result['targets']['M']['transcripts'] == {'ENSTM3': 400}
    assert result['targets']['M']['chosen_by'] == 'Ensembl_canonical'
    assert set(result['targets']['P']['tm']) == {'binding_tm'}


def test_check_flags_probes_off_the_reference_transcript(toy, capsys):
    # transcripts: all also tiles ENSTP2's alternative last exon, which the MANE mRNA lacks.
    table, summary = design(toy, 'all')
    reference = ['--transcripts', str(toy / 'transcript.fa'), '--transcript', 'P=ENSTP1', '--transcript', 'M=ENSTM3']
    assert check_probes.main([str(table), *reference, '--json']) == 1
    result = json.loads(capsys.readouterr().out)
    off = [f for f in result['failures'] if 'not on the designed transcript (ENSTP1)' in f]
    assert off and all(f.startswith('P_') for f in off)
    # With the run's own summary the windows lie on its transcripts, but the choice is flagged.
    assert check_probes.main([str(table), '--summary', str(summary), '--transcripts', str(toy / 'transcript.fa')]) == 0
    assert 'designed on all annotated transcripts' in capsys.readouterr().out
