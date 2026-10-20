"""The U-Probe Agent Skill and its installer.

``u-probe/`` next to this module is an Agent Skill (https://agentskills.io/specification):
``SKILL.md`` with ``name`` and ``description`` frontmatter, ``references/`` and ``scripts/``.
It ships with the package; ``uprobe skill path`` prints where, and ``uprobe skill install``
copies it to where a coding agent discovers skills.

Discovery locations, as documented by each agent:

- Codex: ``.agents/skills/`` in the working directory and its parents up to the repository
  root, and ``$HOME/.agents/skills`` for the user
  (https://developers.openai.com/codex/skills).
- Claude Code: ``.claude/skills/`` in the project (the start directory and its parents up to
  the repository root), and ``~/.claude/skills/`` for the user
  (https://code.claude.com/docs/en/skills).
"""
import shutil
import tempfile
from pathlib import Path

import yaml

SKILL_NAME = 'u-probe'
# Skills directory of each agent, relative to the project directory or to the home directory.
AGENT_SKILL_DIRS = {
    'codex': Path('.agents') / 'skills',
    'claude': Path('.claude') / 'skills',
}
_IGNORE = shutil.ignore_patterns('__pycache__', '*.pyc', '*.pyo', '.DS_Store')


class SkillInstallError(Exception):
    """The skill cannot be installed at the requested place."""


def skill_path() -> Path:
    """The skill folder shipped with this installation of U-Probe."""
    return Path(__file__).resolve().parent / SKILL_NAME


def read_frontmatter(skill_md) -> dict:
    """The YAML frontmatter of a SKILL.md (the block between the first two ``---`` lines)."""
    lines = Path(skill_md).read_text(encoding='utf-8').splitlines()
    if not lines or lines[0].strip() != '---':
        raise ValueError(f"{skill_md}: no YAML frontmatter (the first line must be '---')")
    try:
        end = next(i for i, line in enumerate(lines[1:], start=1) if line.strip() == '---')
    except StopIteration:
        raise ValueError(f"{skill_md}: the frontmatter has no closing '---'") from None
    data = yaml.safe_load('\n'.join(lines[1:end])) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{skill_md}: the frontmatter is not a mapping")
    return data


def skills_dir(target: str, user: bool = False, project_dir=None, home=None) -> Path:
    """The directory that receives the ``u-probe`` folder.

    ``target`` is ``codex``, ``claude`` or a skills directory. For an agent, ``user`` picks
    the user-level directory under ``home`` (default: the home directory) instead of the
    project-level one under ``project_dir`` (default: the working directory).
    """
    if target in AGENT_SKILL_DIRS:
        base = Path(home) if home is not None else Path.home()
        if not user:
            base = Path(project_dir) if project_dir is not None else Path.cwd()
        return base / AGENT_SKILL_DIRS[target]
    if user:
        raise SkillInstallError(f"--user applies to 'codex' and 'claude' only; {target!r} is a directory "
                                "and is used as given.")
    return Path(target).expanduser()


def _is_uprobe_skill(path: Path) -> bool:
    try:
        return read_frontmatter(path / 'SKILL.md').get('name') == SKILL_NAME
    except (OSError, ValueError, yaml.YAMLError):
        return False


def install_skill(directory, force: bool = False) -> Path:
    """Copy the skill into ``directory/u-probe`` and return that path.

    An existing ``u-probe`` folder is replaced only with ``force``, and only if it is a
    U-Probe skill (or a symbolic link, which is removed, not followed).
    """
    source = skill_path()
    if not (source / 'SKILL.md').is_file():
        raise SkillInstallError(f"The skill is missing from this installation ({source}); reinstall U-Probe.")
    directory = Path(directory).expanduser().absolute()
    dest = directory / SKILL_NAME
    if dest.resolve() == source:
        raise SkillInstallError(f"{dest} is the packaged skill itself.")
    replace = dest.exists() or dest.is_symlink()
    if replace:
        if not force:
            raise SkillInstallError(f"{dest} exists; use --force to replace it.")
        if not dest.is_symlink() and not _is_uprobe_skill(dest):
            raise SkillInstallError(f"{dest} exists and is not a U-Probe skill; not replacing it.")
    directory.mkdir(parents=True, exist_ok=True)
    # Copy next to the destination first, so a failed copy leaves any old skill in place.
    staging = Path(tempfile.mkdtemp(prefix=f'.{SKILL_NAME}-', dir=directory))
    try:
        copied = staging / SKILL_NAME
        shutil.copytree(source, copied, ignore=_IGNORE)
        if replace:
            if dest.is_symlink():
                dest.unlink()
            else:
                shutil.rmtree(dest)
        copied.rename(dest)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return dest
