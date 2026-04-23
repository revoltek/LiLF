# Configuration loader for the LiLF pipeline.
# Reads an optional LiLF.conf file from the current or parent directory and
# fills in all pipeline sections with their default values.

import sys
import configparser
import re
import ast
import os
import glob

from LiLF.lib_log import logger


# ---------------------------------------------------------------------------
# Default parameters per kind (and global)
# ---------------------------------------------------------------------------

DEFAULTS = {
    'global': {
        'working_dir': os.getcwd(),
        'max_workers': 1,
        'dry_run': False,
    },
    'download': {
        'download_file': None,
        'macaroon_file': None,
    },
    'setup': {
        'fix_table':              True,
        'keep_IS':                True,
        'backup_full_res':        False,
        'demix_sources':          '',
        'demix_skymodel':         '',
        'demix_field_skymodel':   'LOTSS-DR3',
        'run_aoflagger':          False,
        'validate':               False,
        'tar':                    True,
    },
    'calibrator': {
        'data_dir':               'data-bkp/',
        'skymodel':               '',
        'imaging':                False,
        'fillmissingedges':       True,
        'less_aggressive_flag':   False,
        'develop':                False,
        'use_GNSS':               False,
        'use_shm':                False,
    },
    'timesplit': {
        'data_dir':               'data-bkp/',
        'cal_dir':                '',
        'ngroups':                1,
        'initc':                  0,
        'fillmissingedges':       True,
        'apply_fr':               False,
        'no_aoflagger':           False,
        'ateam_clip':             '',
        'use_GNSS':               False,
    },
    'combine': {},
    'ddparallel': {
        'maxIter':                2,
        'subfield':               '',
        'subfield_min_flux':      20,
        'ph_sol_mode':            'phase',
        'remove3c':               True,
        'fulljones':              False,
        'min_facets':             '',
        'max_facets':             '',
        'develop':                False,
        'data_dir':               '',
        'use_shm':                False,
    },
    'ddserial': {
        'maxIter':                1,
        'minCalFlux60':           0.8,
        'solve_amp':              True,
        'use_shm':                False,
        'use_shm_ddcal':          True,
        'target_dir':             '',
        'manual_dd_cal':          '',
        'develop':                False,
    },
    'quality': {},
    'extract': {
        'max_niter':              10,
        'subtract_region':        '',
        'ph_sol_mode':            'phase',
        'amp_sol_mode':           'diagonal',
        'beam_cut':               0.3,
        'no_selfcal':             False,
        'ampcal':                 'auto',
        'extractRegion':          'target.reg',
    },
    'facetselfcal': {},
}


# ---------------------------------------------------------------------------
# Mandatory parameters per kind
# A step will fail validation if any of these are None or empty string.
# ---------------------------------------------------------------------------

MANDATORY = {
    'global':       [],
    #'download':     ['project', 'obsid', 'beam', 'output'],
    'download':     ['output'],
    'setup':        ['input', 'output'],
    'calibrator':   ['input', 'output'],
    'timesplit':    ['input_h5', 'input_mss', 'output'],
    'combine':      ['input', 'output'],
    'ddparallel':   ['input', 'output'],
    'ddserial':     ['input', 'output'],
    'quality':      ['input'],
    'extract':      ['input', 'output', 'extractRegion'],
    'facetselfcal': ['input', 'output'],
}


def get_kind_defaults(kind: str) -> dict:
    return dict(DEFAULTS.get(kind, {}))


def get_mandatory(kind: str) -> list:
    return list(MANDATORY.get(kind, []))


# ---------------------------------------------------------------------------
# Value parser
# ---------------------------------------------------------------------------

def parse_value(val):
    if val is None or val.strip() == '':
        return None
    v = val.strip()
    if v.lower() in ('true', 'yes'):  return True
    if v.lower() in ('false', 'no'): return False
    if v.startswith('['):
        try:
            return ast.literal_eval(v)
        except Exception:
            inner = v.strip('[]')
            return [x.strip() for x in inner.split(',')]
    try:    return int(v)
    except: pass
    try:    return float(v)
    except: pass
    return v


def quote_list_identifiers(s: str) -> str:
    return re.sub(r"(?<!['\"])(\b[a-zA-Z_]\w*\b)(?!['\"])", r"'\1'", s)


def parse_steps(val: str):
    """
    Accept: [dl1, dl2], cal1, [ts1, ts2] 
    and return a flat-top nested list of strings.
    """
    v = val.strip()
    quoted = quote_list_identifiers(f'[{v}]')
    return ast.literal_eval(quoted)


# ---------------------------------------------------------------------------
# Validation result
# ---------------------------------------------------------------------------

class ValidationError:
    def __init__(self, step_name, kind, missing):
        self.step_name = step_name
        self.kind      = kind
        self.missing   = missing   # list of missing param names

    def __repr__(self):
        return (f"ValidationError(step={self.step_name!r}, "
                f"kind={self.kind!r}, missing={self.missing})")


# ---------------------------------------------------------------------------
# Step
# ---------------------------------------------------------------------------

class Step:
    """
    One pipeline step.

    Parameter resolution order (highest to lowest priority):
        1. local params    (step's own [section] in the parset)
        2. kind defaults   (DEFAULTS[kind])
        3. global params   (top-level parset keys + DEFAULTS['global'])
    """

    def __init__(self, name, kind, local_params, global_params):
        self.name           = name
        self.kind           = kind
        self._local         = local_params
        self._kind_defaults = get_kind_defaults(kind)
        self._mandatory     = get_mandatory(kind)
        self._global        = global_params

        # Auto-set parset_dir unless the user already provided one.
        # Prefer a node-specific parsets/<hostname>/<kind> directory so that
        # different machines can carry custom DP3/losoto parsets without
        # modifying the shared defaults.  Fall back to parsets/<kind>.
        if 'parset_dir' not in self._local or self._local['parset_dir'] is None:
            _lilf_dir  = os.path.dirname(os.path.dirname(__file__))
            _parset_dir = os.path.join(_lilf_dir, 'parsets', kind)
            if os.path.isdir(_parset_dir):
                self._local['parset_dir'] = _parset_dir

    # --- param access -------------------------------------------------------

    def __getitem__(self, key):
        if key in self._local and self._local[key] is not None:
            return self._local[key]
        if key in self._kind_defaults and self._kind_defaults[key] is not None:
            return self._kind_defaults[key]
        return self._global.get(key)

    def get(self, key, default=None):
        val = self[key]
        return val if val is not None else default

    def all_params(self):
        merged = {**self._global}
        merged.update({k: v for k, v in self._kind_defaults.items()
                        if v is not None})
        merged.update({k: v for k, v in self._local.items()
                        if v is not None})
        return merged

    # --- validation ---------------------------------------------------------

    def validate(self):
        """
        Check that all mandatory parameters for this step's kind are
        set and non-empty.  Returns a ValidationError or None.
        """
        missing = []
        for key in self._mandatory:
            val = self[key]
            if val is None or val == '':
                missing.append(key)
        if missing:
            return ValidationError(self.name, self.kind, missing)
        return None

    def is_valid(self):
        return self.validate() is None

    # --- display ------------------------------------------------------------

    def describe(self):
        lines = [f'  Step: {self.name}  (kind={self.kind})']

        all_keys = (set(self._global)
                    | set(self._kind_defaults)
                    | set(self._local))

        for k in sorted(all_keys):
            if k in self._local and self._local[k] is not None:
                src = 'local'
                v   = self._local[k]
            elif k in self._kind_defaults and self._kind_defaults[k] is not None:
                src = 'default'
                v   = self._kind_defaults[k]
            else:
                v   = self._global.get(k)
                src = 'global'
            if v is not None:
                tag = ' *' if k in self._mandatory else '  '
                lines.append(f'  {tag} {k:30s} = {str(v):30s}  [{src}]')

        err = self.validate()
        if err:
            lines.append(f'  !! MISSING MANDATORY: {err.missing}')

        return '\n'.join(lines)

    def __repr__(self):
        return (f"Step(name={self.name!r}, kind={self.kind!r}, "
                f"params={self.all_params()})")


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

class Pipeline:
    """
    Parsed pipeline.

    execution_plan() returns an ordered list where each entry is either:
        Step          -> run serially, block until done
        list[Step]    -> run in parallel, block until ALL done
    """

    def __init__(self, global_params, steps_order, step_objects):
        self.global_params = global_params
        self.steps_order   = steps_order
        self.steps         = step_objects

    # --- plan ---------------------------------------------------------------

    def execution_plan(self):
        plan = []
        for item in self.steps_order:
            if isinstance(item, list):
                group = [self.steps[n] for n in item if n in self.steps]
                if group:
                    plan.append(group)
            else:
                if item in self.steps:
                    plan.append(self.steps[item])
        return plan

    def flat_steps(self):
        result = []
        for stage in self.execution_plan():
            if isinstance(stage, list):
                result.extend(stage)
            else:
                result.append(stage)
        return result

    # --- queries ------------------------------------------------------------

    def steps_of_kind(self, kind):
        return [s for s in self.steps.values() if s.kind == kind]

    def parallel_groups(self):
        return [s for s in self.execution_plan() if isinstance(s, list)]

    def serial_steps(self):
        return [s for s in self.execution_plan() if isinstance(s, Step)]

    # --- validation ---------------------------------------------------------

    def validate(self):
        """
        Validate all steps. Returns a list of ValidationError (empty = all ok).
        """
        errors = []
        for step in self.flat_steps():
            err = step.validate()
            if err:
                errors.append(err)
        return errors

    def validate_and_report(self):
        errors = self.validate()
        if not errors:
            logger.info('Validation OK — all mandatory parameters are set.')
        else:
            logger.error(f'Validation FAILED — {len(errors)} step(s) have missing mandatory params:')
            for e in errors:
                logger.error(f'  [{e.step_name}] ({e.kind}) missing: {e.missing}')
        return errors

    def get(self, key, default=None):
        return self.global_params.get(key, default)

    # --- display ------------------------------------------------------------

    def describe(self, verbose=False):
        lines = ['\n=== Global Parameters ===']
        for k, v in self.global_params.items():
            lines.append(f'  {k} = {v}')

        lines.append('\n=== Execution Plan ===')
        for i, stage in enumerate(self.execution_plan(), 1):
            if isinstance(stage, list):
                entries = ', '.join(f'{s.name}({s.kind})' for s in stage)
                valid   = all(s.is_valid() for s in stage)
                flag    = '' if valid else '  !! INVALID'
                lines.append(f'  [{i:02d}] PARALLEL  {{ {entries} }}{flag}')
                if verbose:
                    for s in stage:
                        lines.append(s.describe())
            else:
                flag = '' if stage.is_valid() else '  !! INVALID'
                lines.append(
                    f'  [{i:02d}] SERIAL    {stage.name}  '
                    f'(kind={stage.kind}){flag}')
                if verbose:
                    lines.append(stage.describe())

        #errors = self.validate()
        #if errors:
        #    lines.append('\n=== Validation Errors ===')
        #    for e in errors:
        #        lines.append(f'  [{e.step_name}] ({e.kind}) missing: {e.missing}')

        return '\n'.join(lines)

    def __repr__(self):
        return self.describe()


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def read_parset(path: str) -> Pipeline:

    if path is None:
        # Auto-discovery: try each conventional extension in order.
        _candidates = [
            f for ext in ('conf', 'config', 'cfg', 'parset')
            for f in glob.glob(f'[Ll][Ii][Ll][Ff].{ext}')
        ]
        if len(_candidates) > 1:
            raise LookupError(f'Found more than one configuration file: {_candidates}')
        elif len(_candidates) == 1:
            path = _candidates[0]
        else:
            logger.error(
            'No config file given and none found in the current directory. '
            'Pass a config file as argument or create LiLF.conf here.'
            )
            sys.exit(1)

    if not os.path.exists(path):
        logger.error(f'Config file not found: {path}')
        sys.exit(1)

    with open(path) as f:
        raw = f.read()

    steps_order = []
    steps_match = re.search(r'^steps\s*=\s*(.+)$', raw, re.MULTILINE)
    if steps_match:
        steps_order = parse_steps(steps_match.group(1))
        raw = raw[:steps_match.start()] + raw[steps_match.end():]

    # Split into the global preamble (before any [section]) and the rest.
    section_pat   = re.compile(r'^\[([^\]]+)\]', re.MULTILINE)
    first_section = section_pat.search(raw)
    global_text   = raw[:first_section.start()] if first_section else raw
    sections_text = raw[first_section.start():]  if first_section else ''

    # Parse global key=value pairs directly — no fake [__global__] section needed.
    global_params = dict(DEFAULTS['global'])
    for m in re.finditer(r'^\s*([A-Za-z_]\w*)\s*=\s*(.*)$', global_text, re.MULTILINE):
        global_params[m.group(1)] = parse_value(m.group(2))

    # Parse named [section] blocks with ConfigParser.
    cfg = configparser.RawConfigParser(allow_no_value=True)
    cfg.optionxform = str
    if sections_text:
        cfg.read_string(sections_text)

    step_objects = {}
    for section in cfg.sections():
        local = {}
        for key, val in cfg.items(section):
            local[key] = parse_value(val)

        kind = local.pop('type', section)
        step_objects[section] = Step(
            name=section,
            kind=kind,
            local_params=local,
            global_params=global_params,
        )

    return Pipeline(global_params, steps_order, step_objects)