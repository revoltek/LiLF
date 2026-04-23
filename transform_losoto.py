#!/usr/bin/env python3
"""Transform run_losoto calls from old to new signature.

Old: lib_scheduler.run_losoto(
        s,
        [h5s],
        parsets,
        logname='losoto-' + (c) + '.log',
        h5_out=os.path.join(None, 'cal-' + (c) + '.h5'),
        plots_dir=None)
New: lib_scheduler.run_losoto(
        s,
        [parsets],
        logname='losoto.log',
        logname='losoto-' + (h5s) + '.log',
        h5_out='cal-' + (h5s) + '.h5',
        plots_dir=None,
        h5_out=None)
"""

import os
import glob


def find_matching_paren(text, start):
    depth = 0
    i = start
    in_str = False
    str_char = None

    while i < len(text):
        ch = text[i]
        if in_str:
            if ch == '\\':
                i += 2
                continue
            if isinstance(str_char, str) and len(str_char) == 3:
                if text[i:i+3] == str_char:
                    in_str = False
                    str_char = None
                    i += 3
                    continue
            elif ch == str_char:
                in_str = False
                str_char = None
        else:
            if ch in ('"', "'"):
                if text[i:i+3] in ('"""', "'''"):
                    str_char = text[i:i+3]
                    in_str = True
                    i += 3
                    continue
                else:
                    str_char = ch
                    in_str = True
            elif ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
                if depth == 0:
                    return i
        i += 1
    return -1


def split_args(s):
    args = []
    depth = 0
    in_str = False
    str_char = None
    current = []
    i = 0

    while i < len(s):
        ch = s[i]
        if in_str:
            current.append(ch)
            if ch == '\\':
                i += 1
                if i < len(s):
                    current.append(s[i])
            elif isinstance(str_char, str) and len(str_char) == 3:
                if s[i:i+3] == str_char:
                    current.append(s[i+1])
                    current.append(s[i+2])
                    i += 2
                    in_str = False
                    str_char = None
            elif ch == str_char:
                in_str = False
                str_char = None
        else:
            if ch in ('"', "'"):
                if s[i:i+3] in ('"""', "'''"):
                    str_char = s[i:i+3]
                    current.append(ch)
                    current.append(s[i+1])
                    current.append(s[i+2])
                    i += 2
                    in_str = True
                else:
                    str_char = ch
                    current.append(ch)
                    in_str = True
            elif ch in ('(', '[', '{'):
                depth += 1
                current.append(ch)
            elif ch in (')', ']', '}'):
                depth -= 1
                current.append(ch)
            elif ch == ',' and depth == 0:
                args.append(''.join(current))
                current = []
                i += 1
                continue
            else:
                current.append(ch)
        i += 1

    if current:
        args.append(''.join(current))
    return args


def derive_logname_h5out(c_expr):
    c = c_expr.strip()

    # Simple single-quoted string: 'xxx'
    if len(c) >= 2 and c[0] == "'" and c[-1] == "'" and c.count("'") == 2:
        inner = c[1:-1]
        if inner.endswith('.h5'):
            stem = os.path.splitext(os.path.basename(inner))[0]
            return f"'losoto-{stem}.log'", c
        else:
            return f"'losoto-{inner}.log'", f"'cal-{inner}.h5'"

    # Simple double-quoted string: "xxx"
    if len(c) >= 2 and c[0] == '"' and c[-1] == '"' and c.count('"') == 2:
        inner = c[1:-1]
        if inner.endswith('.h5'):
            stem = os.path.splitext(os.path.basename(inner))[0]
            return f"'losoto-{stem}.log'", c
        else:
            return f"'losoto-{inner}.log'", f"'cal-{inner}.h5'"

    # f-string: f'xxx' or f"xxx"
    if (c.startswith("f'") and c.endswith("'")) or (c.startswith('f"') and c.endswith('"')):
        inner = c[2:-1]
        if inner.endswith('.h5'):
            last_slash = inner.rfind('/')
            last_seg = inner[last_slash+1:] if last_slash >= 0 else inner
            stem = last_seg[:-3]
            return f"'losoto-{stem}.log'", c
        else:
            q = c[1]  # quote char
            return f"f{q}losoto-{inner}.log{q}", f"f{q}cal-{inner}.h5{q}"

    # Generic expression (string concat, % format, variable, etc.)
    # Assume does NOT end in .h5
    return f"'losoto-' + ({c}) + '.log'", f"'cal-' + ({c}) + '.h5'"


def transform_call(full_call, ws_indent):
    paren_start = full_call.index('(')
    inner = full_call[paren_start+1:-1]

    args = split_args(inner)

    if len(args) < 4:
        return None  # can't transform

    arg_s = args[0].strip()
    arg_c = args[1].strip()
    arg_h5s = args[2].strip()
    arg_parsets = args[3].strip()
    rest = [a.strip() for a in args[4:]]

    # Remove trailing backslash-newlines from args (line continuation inside parens)
    def clean_arg(a):
        return a.replace('\\\n', '').strip()

    arg_c = clean_arg(arg_c)
    arg_h5s = clean_arg(arg_h5s)
    arg_parsets = clean_arg(arg_parsets)
    rest = [clean_arg(a) for a in rest if clean_arg(a)]

    # Separate h5_dir
    h5_dir_val = None
    other_kwargs = []
    for kw in rest:
        if kw.startswith('h5_dir='):
            h5_dir_val = kw[7:].strip()
        else:
            other_kwargs.append(kw)

    logname, h5_out = derive_logname_h5out(arg_c)

    if h5_dir_val:
        if h5_out.startswith("'") and h5_out.endswith("'"):
            bn = os.path.basename(h5_out[1:-1])
            h5_out = f"os.path.join({h5_dir_val}, '{bn}')"
        else:
            h5_out = f"os.path.join({h5_dir_val}, os.path.basename({h5_out}))"

    # Wrap h5s in list if not already a list
    if not arg_h5s.startswith('['):
        arg_h5s = f'[{arg_h5s}]'

    new_args = [arg_s, arg_h5s, arg_parsets,
                f'logname={logname}', f'h5_out={h5_out}'] + other_kwargs

    # Format: use multi-line
    cont_indent = ws_indent + '        '
    args_str = f',\n{cont_indent}'.join(new_args)
    return f'lib_scheduler.run_losoto(\n{cont_indent}{args_str})'


def transform_file(content):
    PATTERN = 'lib_scheduler.run_losoto('
    result = []
    i = 0

    while i < len(content):
        idx = content.find(PATTERN, i)
        if idx == -1:
            result.append(content[i:])
            break

        # Check if in a comment
        line_start = content.rfind('\n', 0, idx) + 1
        line_prefix = content[line_start:idx]
        is_comment = line_prefix.lstrip().startswith('#')

        result.append(content[i:idx])

        paren_pos = idx + len(PATTERN) - 1
        close_paren = find_matching_paren(content, paren_pos)

        if close_paren == -1:
            result.append(PATTERN)
            i = idx + len(PATTERN)
            continue

        full_call = content[idx:close_paren+1]

        if is_comment:
            result.append(full_call)
            i = close_paren + 1
            continue

        # Get whitespace indent
        ws_indent = ''
        for ch in line_prefix:
            if ch in (' ', '\t'):
                ws_indent += ch
            else:
                break

        new_call = transform_call(full_call, ws_indent)
        result.append(new_call if new_call is not None else full_call)
        i = close_paren + 1

    return ''.join(result)


base = '/home/baq1889/LiLF'
py_files = sorted(glob.glob(f'{base}/**/*.py', recursive=True))

changed = []
errors = []
for filepath in py_files:
    with open(filepath, 'r') as f:
        content = f.read()

    if 'lib_scheduler.run_losoto(' not in content:
        continue

    try:
        new_content = transform_file(content)
        if new_content != content:
            with open(filepath, 'w') as f:
                f.write(new_content)
            changed.append(os.path.relpath(filepath, base))
            print(f'Modified: {os.path.relpath(filepath, base)}')
    except Exception as e:
        errors.append((filepath, str(e)))
        print(f'ERROR in {filepath}: {e}')

print(f'\nModified {len(changed)} files')
if errors:
    print(f'Errors: {len(errors)}')
    for f, e in errors:
        print(f'  {f}: {e}')
