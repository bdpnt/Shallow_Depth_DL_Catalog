# PATTERNS.md — Python Coding Conventions

This file documents the coding conventions used throughout the Shallow_Depth_DL_Catalog project. Follow these patterns when writing new Python scripts.

---

## Module docstring

Every script starts with a docstring in this format:

```
"""
name.py
============================
One-line purpose.

Longer explanation of what the script does.

Usage
-----
    python package/name.py --arg value

    # Minimal / all defaults
    python package/name.py
"""
```

---

## Import order

1. stdlib (alphabetical)
2. blank line
3. third-party (alphabetical)
4. blank line
5. local — when importing from another sub-package, add a `sys.path.insert` guard before the import so the script works both when run directly and when imported as part of the package:

```python
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from subpackage.module import func
```

---

## Module-level path constants

```python
_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)   # adjust depth as needed
```

Default input/output paths are built from these anchors with `os.path.join`. The exact paths depend on what the script consumes and produces — there is no fixed convention for the values, only for how they are derived: always from `__file__`, never hardcoded strings.

---

## Logger setup

Scripts that produce persistent diagnostic output use the `logging` module with a `FileHandler` only — no console output.

```python
logger = logging.getLogger('module_name')   # module-level singleton

def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)
    return log_path
```

Pure utility or visualization scripts that do not need persistent logs use plain `print()`.

---

## Section separators

Use section separators to group related code whenever a script has distinct logical blocks:

```python
# ---------------------------------------------------------------------------
# Section name
# ---------------------------------------------------------------------------
```

Two sections are **mandatory** in every script:
- `# Public API` — contains the callable function(s) exposed to other modules
- `# CLI entry point` — contains `main()` and the `if __name__` guard

Other sections (e.g. `# Internal helpers`, `# Constants`, `# Format handlers`) are used as needed when the script is large enough to benefit from them.

---

## Public API + CLI pattern

Every script exposes:
1. A named public function — the real implementation.
2. A thin `main()` that parses `argparse` and calls it.
3. `if __name__ == '__main__': main()` at the bottom.

```python
def do_thing(input_path, output_path=None, log_dir=None):
    """NumPy-style docstring."""
    ...
    return {'output': output_path, 'log': log_path, ...}


def main():
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('--input',  required=True)
    parser.add_argument('--output', default=None)
    args = parser.parse_args()
    do_thing(args.input, args.output)


if __name__ == '__main__':
    main()
```

This makes every script importable (`from package.module import do_thing`) and runnable from the command line without any change.

---

## Return value convention

Public API functions always return a summary `dict`. Mandatory keys:

- `'output'` — path to the output file
- `'log'` — path to the log file (when logging is used)
- Counts as needed: `n_events`, `n_picks`, `n_added`, `n_skipped_*`, …

---

## Docstrings

NumPy style, with `Parameters`, `Returns`, and `Notes` sections as needed.

---

## Private helpers

Internal helpers are prefixed with `_` (e.g. `_parse_event_header`, `_setup_logger`).

---

## Extensibility pattern

Use a dispatch dict for format-specific or variant-specific logic:

```python
FORMAT_HANDLERS = {
    'FORMAT_A': _handle_format_a,
    # Add new formats here as handler functions
}
```

---

## Data classes

Use `@dataclass` for structured data containers. List fields use `field(default_factory=list)`.

---

## Figures

Use `seaborn` to style figures — apply a theme and palette at the top of the plotting function. Always call `plt.close(fig)` after saving to avoid memory leaks.

---

## Scope when applying these conventions

These conventions cover **structure and style only**: module layout, import order, path constants, logger setup, section separators, public API shape, and return dicts.

They do **not** authorise silently changing:
- Data transformation logic (column indices, filters, arithmetic)
- Algorithm or matching logic
- Any existing behaviour that is already correct

When applying PATTERNS.md to an existing script, reformat and restructure. If something in the existing logic looks wrong or improvable, **ask the user before changing it** — do not alter it silently. The user can then verify and confirm or reject the change.
