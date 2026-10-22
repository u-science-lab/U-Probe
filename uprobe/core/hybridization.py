"""Hybridization conditions for Tm (primer3) and fold_score (ViennaRNA).

A protocol may declare the conditions its probes hybridize under::

    hybridization:
      temperature_c: 37    # hybridization temperature (°C), used by fold_score
      na_mM: 390           # monovalent cations, e.g. 2x SSC
      mg_mM: 0             # divalent cations
      dntp_mM: 0           # dNTPs (chelate Mg2+; no effect when mg_mM is 0)
      probe_nM: 50         # concentration of each oligo
      formamide_pct: 30    # formamide, % (v/v)

An ``annealing_temperature`` or ``fold_score`` attribute may carry its own
``hybridization`` mapping whose fields override the protocol-level ones.
Fields left unset keep the legacy values (the fisheye-compatible primer3
arguments and ViennaRNA's 37 °C), so protocols without the block produce
unchanged results.

An ``annealing_temperature`` attribute computes Tm with primer3 unless it sets
``tm_method`` to the formula of a published design code, under the same conditions:
``multi_padlock_design`` (HybISS padlocks, cal_temp_multi_padlock_design),
``merfish_design`` (MERFISH, cal_temp_merfish_design) or ``oligominer``
(OligoMiner, cal_temp_oligominer).
"""
import math


FIELDS = ('temperature_c', 'na_mM', 'mg_mM', 'dntp_mM', 'probe_nM', 'formamide_pct')
ATTRIBUTE_TYPES = ('annealing_temperature', 'fold_score')
# Legacy conditions: the defaults of cal_temp and of ViennaRNA's model.
TM_DEFAULTS = {'na_mM': 50.0, 'mg_mM': 0.0, 'dntp_mM': 0.6, 'probe_nM': 50.0, 'formamide_pct': 0.0}
FOLD_DEFAULTS = {'temperature_c': 37.0}
# Tm formulas of annealing_temperature attributes (option tm_method); the first is the default.
TM_METHODS = ('primer3', 'multi_padlock_design', 'merfish_design', 'oligominer')
# Hybridization fields each design-code formula uses (the others are ignored); primer3 uses all but temperature_c.
TM_FORMULA_FIELDS = {'multi_padlock_design': ('na_mM', 'mg_mM', 'probe_nM', 'formamide_pct'),
                     'merfish_design': ('na_mM', 'probe_nM'),
                     'oligominer': ('na_mM', 'probe_nM', 'formamide_pct')}
PRIMER3_ARGUMENTS = {'na_mM': 'mv_conc', 'mg_mM': 'dv_conc', 'dntp_mM': 'dntp_conc', 'probe_nM': 'dna_conc'}
# primer3 lowers Tm by (2.88 - 0.453 * GC fraction) °C per mol/l formamide
# (Blake & Delcourt 1996, Nucleic Acids Res 24:2095-2103,
# doi:10.1093/nar/24.11.2095).  Protocols give % (v/v), converted with
# formamide's density (1.1334 g/ml at 20 °C) and molar mass (45.041 g/mol),
# both from PubChem CID 713: 1% = 0.2516 mol/l, i.e. 0.61-0.72 °C per %,
# in line with the linear 0.60 °C/% of Hutton 1977 (Nucleic Acids Res
# 4:3537-3555, doi:10.1093/nar/4.10.3537).
FORMAMIDE_MOLAR_PER_PCT = 10 * 1.1334 / 45.041
# Auto-generated Tm filters: 50-90 °C under the legacy conditions; under
# declared conditions [T, T + 40] °C at the hybridization temperature T.
LEGACY_TM_WINDOW = (50.0, 90.0)
TM_WINDOW_WIDTH = LEGACY_TM_WINDOW[1] - LEGACY_TM_WINDOW[0]


def validate_hybridization(block, where='hybridization'):
    """Return the fields of one hybridization mapping as floats."""
    if block is None:
        return {}
    if not isinstance(block, dict):
        raise ValueError(f'{where} must be a mapping with fields: {", ".join(FIELDS)}')
    unknown = [str(key) for key in block if key not in FIELDS]
    if unknown:
        raise ValueError(f'{where}: unknown field(s) {", ".join(unknown)}; expected: {", ".join(FIELDS)}')
    values = {}
    for key, value in block.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f'{where}.{key} must be a number, got {value!r}')
        if value < 0:
            raise ValueError(f'{where}.{key} must be >= 0, got {value}')
        values[key] = float(value)
    for key in ('na_mM', 'probe_nM'):
        if values.get(key) == 0:
            raise ValueError(f'{where}.{key} must be > 0')
    for key in ('temperature_c', 'formamide_pct'):
        if values.get(key, 0) > 100:
            raise ValueError(f'{where}.{key} must be between 0 and 100, got {values[key]:g}')
    return values


def resolve_conditions(protocol, attr_name):
    """Declared conditions for one attribute: protocol block, then its override."""
    conditions = validate_hybridization(protocol.get('hybridization'))
    attr = (protocol.get('attributes') or {}).get(attr_name)
    if isinstance(attr, dict):
        conditions.update(validate_hybridization(
            attr.get('hybridization'), f'attributes.{attr_name}.hybridization'))
    return conditions


def tm_arguments(conditions):
    """cal_temp keyword arguments for the declared Tm fields only."""
    arguments = {PRIMER3_ARGUMENTS[key]: value for key, value in conditions.items()
                 if key in PRIMER3_ARGUMENTS}
    if 'formamide_pct' in conditions:
        arguments['formamide_conc'] = conditions['formamide_pct'] * FORMAMIDE_MOLAR_PER_PCT
    return arguments


def validate_protocol_hybridization(protocol):
    """Check the protocol block and every attribute-level override."""
    validate_hybridization(protocol.get('hybridization'))
    for name, attr in (protocol.get('attributes') or {}).items():
        if not isinstance(attr, dict) or 'hybridization' not in attr:
            continue
        if str(attr.get('type', '')).lower() not in ATTRIBUTE_TYPES:
            raise ValueError(f'attributes.{name}.hybridization is only supported for '
                             f'{" and ".join(ATTRIBUTE_TYPES)} attributes')
        resolve_conditions(protocol, name)
    for name, attr in (protocol.get('attributes') or {}).items():
        if isinstance(attr, dict) and 'tm_method' in attr:
            tm_method(attr, f'attributes.{name}')


def tm_method(attr, where='attribute'):
    """The Tm formula of an annealing_temperature attribute: primer3 (default) or a design code's."""
    if not isinstance(attr, dict) or 'tm_method' not in attr:
        return TM_METHODS[0]
    if str(attr.get('type', '')).lower() != 'annealing_temperature':
        raise ValueError(f'{where}.tm_method is only supported for annealing_temperature attributes')
    method = attr['tm_method']
    if method not in TM_METHODS:
        raise ValueError(f'{where}.tm_method must be one of {", ".join(TM_METHODS)}, got {method!r}')
    return method


def tm_formula_arguments(method, conditions):
    """Keyword arguments of the formula ``method``: declared conditions, legacy defaults for the rest."""
    return {key: conditions.get(key, TM_DEFAULTS[key]) for key in TM_FORMULA_FIELDS[method]}


def unused_fields(method, conditions):
    """Declared Tm fields (not temperature_c) that the formula of ``method`` does not use.

    Not reported for multi_padlock_design, whose protocols have always declared dntp_mM.
    """
    if method not in TM_FORMULA_FIELDS or method == 'multi_padlock_design':
        return []
    return [key for key in conditions if key in TM_DEFAULTS and key not in TM_FORMULA_FIELDS[method]]


def multi_padlock_design_arguments(conditions):
    """cal_temp_multi_padlock_design keyword arguments: declared conditions, legacy defaults for the rest.

    The formula has no dNTP term, so dntp_mM is not used.
    """
    return {key: conditions.get(key, TM_DEFAULTS[key]) for key in ('na_mM', 'mg_mM', 'probe_nM', 'formamide_pct')}


def declares_hybridization(protocol):
    return 'hybridization' in protocol or any(
        isinstance(attr, dict) and 'hybridization' in attr
        for attr in (protocol.get('attributes') or {}).values())


def effective_conditions(protocol):
    """Conditions in effect for each Tm/fold attribute, legacy defaults filled in."""
    result = {}
    for name, attr in (protocol.get('attributes') or {}).items():
        kind = str(attr.get('type', '')).lower() if isinstance(attr, dict) else ''
        if kind in ATTRIBUTE_TYPES:
            defaults = TM_DEFAULTS if kind == 'annealing_temperature' else FOLD_DEFAULTS
            declared = resolve_conditions(protocol, name)
            result[name] = {key: declared.get(key, value) for key, value in defaults.items()}
    return result


def default_tm_window(protocol, attr_name):
    """(low, high, temperature source) of the auto-generated Tm filter of one attribute.

    The window is [T, T + 40] °C at the hybridization temperature T of the
    conditions the attribute's Tm is computed under: its own temperature_c
    ('attribute'), else the protocol's ('protocol'). An override that changes
    the salt, probe or formamide conditions without its own temperature_c
    (e.g. the design conditions of a published Tm window) is not under the
    protocol's conditions and keeps the legacy 50-90 °C (source None), as
    does a protocol without declared Tm conditions. Tm conditions declared
    without any temperature use the default 37 °C ('default').
    """
    block = validate_hybridization(protocol.get('hybridization'))
    attr = (protocol.get('attributes') or {}).get(attr_name)
    override = validate_hybridization(attr.get('hybridization'), f'attributes.{attr_name}.hybridization') \
        if isinstance(attr, dict) else {}
    own_conditions = any(override[key] != block.get(key, default)
                         for key, default in TM_DEFAULTS.items() if key in override)
    if 'temperature_c' in override:
        temperature, source = override['temperature_c'], 'attribute'
    elif own_conditions:
        return LEGACY_TM_WINDOW + (None,)
    elif 'temperature_c' in block:
        temperature, source = block['temperature_c'], 'protocol'
    elif any(key in block for key in TM_DEFAULTS):
        temperature, source = FOLD_DEFAULTS['temperature_c'], 'default'
    else:
        return LEGACY_TM_WINDOW + (None,)
    return temperature, temperature + TM_WINDOW_WIDTH, source


def describe_conditions(protocol):
    """One line per distinct condition set, listing the attributes it applies to."""
    units = {'temperature_c': '{:g} °C', 'na_mM': '{:g} mM Na+', 'mg_mM': '{:g} mM Mg2+',
             'dntp_mM': '{:g} mM dNTP', 'probe_nM': '{:g} nM probe', 'formamide_pct': '{:g}% formamide'}
    groups = {}
    attributes = protocol.get('attributes') or {}
    for name, conditions in effective_conditions(protocol).items():
        kind = 'Tm' if 'na_mM' in conditions else 'fold'
        method = tm_method(attributes.get(name), f'attributes.{name}') if kind == 'Tm' else TM_METHODS[0]
        if method != TM_METHODS[0]:
            kind = f'Tm [{method}]'
            conditions = {key: value for key, value in conditions.items() if key in TM_FORMULA_FIELDS[method]}
        text = ', '.join(units[key].format(value) for key, value in conditions.items())
        groups.setdefault((kind, text), []).append(name)
    return [f'{kind} ({", ".join(names)}): {text}' for (kind, text), names in groups.items()]
