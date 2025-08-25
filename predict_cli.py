#!/usr/bin/env python3
"""
Command-line predictor that shares the same enhanced recommendation engine and export builder as the GUI.

Usage (PowerShell example):
  python predict_cli.py '{"reaction_smiles":"Brc1cc...","selected_reaction_type":"C-N Coupling - Ullmann"}'

Or read from stdin:
  Get-Content input.json | python predict_cli.py
"""
from __future__ import annotations

import sys
import os
import json

# Ensure project root on path
_HERE = os.path.dirname(__file__)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from enhanced_recommendation_engine import create_recommendation_engine
from prediction_export import build_export_payload


def _write_json(obj: dict) -> None:
    """Write JSON to stdout with UTF-8 safety on Windows consoles.

    Tries text write first; on UnicodeEncodeError falls back to writing UTF-8 bytes.
    """
    try:
        # Prefer UTF-8 where possible
        try:
            # Python 3.7+: reconfigure stream to utf-8 if supported
            sys.stdout.reconfigure(encoding='utf-8')  # type: ignore[attr-defined]
        except Exception:
            pass
        s = json.dumps(obj, ensure_ascii=False)
        try:
            sys.stdout.write(s)
        except UnicodeEncodeError:
            # Fallback: write UTF-8 bytes directly
            sys.stdout.buffer.write(s.encode('utf-8'))
    except Exception:
        # Last-resort ASCII-safe output
        sys.stdout.write(json.dumps(obj, ensure_ascii=True))


def _load_input() -> dict:
    """Load input JSON from one of:
    - --input-file/-f <path>
    - JSON passed as CLI args (joined)
    - stdin (only if not a TTY)
    """
    args = sys.argv[1:]

    # Parse known flags (-f/--input-file, -o/--output-file, QUARC flags) and filter others for JSON
    file_path: str | None = None
    filtered: list[str] = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in ("-f", "--input-file") and i + 1 < len(args):
            file_path = args[i + 1]
            i += 2
            continue
        if a in ("-o", "--output-file") and i + 1 < len(args):
            # skip output flag from JSON join
            i += 2
            continue
        if a in ("--use-quarc", "--quarc-config", "--quarc-topk"):
            # skip quarc flags and their values
            i += 2 if (i + 1) < len(args) else 1
            continue
        filtered.append(a)
        i += 1

    # File flag support
    if file_path:
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                return json.loads(f.read())
        except Exception:
            pass

    # Direct single-arg JSON
    if filtered and len(filtered) == 1 and filtered[0].lstrip().startswith('{'):
        try:
            return json.loads(filtered[0])
        except Exception:
            pass

    # Join remaining non-flag args and try parse as JSON
    if filtered:
        try:
            joined = " ".join(filtered).strip()
            if joined:
                return json.loads(joined)
        except Exception:
            # fall through to stdin
            pass

    # Extract JSON between outermost braces from the raw argv string (handles tokenized JSON)
    try:
        raw = " ".join(args)
        start = raw.find('{')
        end = raw.rfind('}')
        if start != -1 and end != -1 and end > start:
            candidate = raw[start:end+1]
            try:
                return json.loads(candidate)
            except Exception:
                # Heuristic: coerce loose dict style into valid JSON
                import re
                loose = candidate
                # Quote keys: key: -> "key":
                loose = re.sub(r'([A-Za-z0-9_\-]+)\s*:', r'"\1":', loose)
                # Quote string values that are unquoted (until comma or closing brace)
                def _qval(m):
                    val = m.group(1).strip()
                    # already quoted or a number/boolean/null
                    if val.startswith('"') or re.fullmatch(r'-?\d+(?:\.\d+)?', val) or val in ('true','false','null'):
                        return ':' + m.group(1)
                    # wrap with quotes, escape existing quotes if any
                    return ': "' + val.replace('"', '\\"') + '"'
                loose = re.sub(r':\s*([^\",}][^,}]*)', _qval, loose)
                return json.loads(loose)
    except Exception:
        pass

    # Environment variable fallback
    try:
        env_payload = os.environ.get('PREDICT_JSON')
        if env_payload:
            return json.loads(env_payload)
    except Exception:
        pass

    # Read stdin only if piped (avoid blocking when interactive)
    try:
        if not sys.stdin.isatty():
            data = sys.stdin.read()
            if data and data.strip():
                return json.loads(data)
    except Exception:
        pass
    return {}


def _extract_output_path() -> str | None:
    """Extract --output-file/-o path from argv if provided."""
    args = sys.argv[1:]
    for i, a in enumerate(args):
        if a in ("-o", "--output-file") and i + 1 < len(args):
            return args[i + 1]
    return None


def _extract_quarc_flags() -> dict:
    """Extract QUARC-related flags from argv to forward into engine options.

    Supported:
      --use-quarc [auto|always|off]
      --quarc-config <path>
      --quarc-topk <N>
    """
    args = sys.argv[1:]
    opts: dict = {}
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--use-quarc" and i + 1 < len(args):
            opts['use_quarc'] = args[i + 1]
            i += 2
            continue
        if a == "--quarc-config" and i + 1 < len(args):
            opts['quarc_config'] = args[i + 1]
            i += 2
            continue
        if a == "--quarc-topk" and i + 1 < len(args):
            try:
                opts['quarc_topk'] = int(args[i + 1])
            except Exception:
                opts['quarc_topk'] = 5
            i += 2
            continue
        i += 1
    return opts


def main() -> int:
    payload = _load_input()
    reaction_smiles = payload.get('reaction_smiles') or ''
    # Normalize reaction type: treat missing/empty/whitespace as Auto-detect
    selected_type = (payload.get('selected_reaction_type') or '').strip() or 'Auto-detect'

    # Guard: no input provided
    if not reaction_smiles:
        _write_json({
            'analysis_type': 'error',
            'status': 'failed',
            'error': 'No input provided. Pass JSON via argument, --input-file, or stdin.'
        })
        return 2

    engine = create_recommendation_engine()
    # If engine supports options, set them via attribute to avoid breaking signature
    try:
        quarc_opts = _extract_quarc_flags()
        if quarc_opts:
            setattr(engine, '_cli_quarc_options', quarc_opts)
    except Exception:
        pass
    recs = engine.get_recommendations(reaction_smiles, selected_type)

    # Normalize to GUI-like result to reuse export builder
    result = {
        'analysis_type': recs.get('analysis_type'),
        'reaction_smiles': reaction_smiles,
        'reaction_type': selected_type,
        'status': recs.get('status', 'success' if 'error' not in recs else 'failed'),
        'recommendations': {
            'reaction_type': recs.get('reaction_type'),
            'ligand_recommendations': recs.get('ligand_recommendations', []),
            'solvent_recommendations': recs.get('solvent_recommendations', []),
            'base_recommendations': recs.get('base_recommendations', []),
            'combined_conditions': recs.get('combined_conditions', []),
            'property_based_alternatives': recs.get('property_based_alternatives', {}),
            'dataset_info': recs.get('dataset_info', {}),
            'related_reactions': recs.get('related_reactions', []),
            'general_recommendations': recs.get('general_recommendations', {}),
        }
    }

    # Also attach general_recommendations at the top level for exporters that look there
    if recs.get('general_recommendations'):
        result['general_recommendations'] = recs['general_recommendations']
    # Propagate providers (e.g., ['quarc_oss','analytics']) if engine set them
    if recs.get('providers'):
        result['providers'] = recs.get('providers')

    export = build_export_payload(result)
    # Attach providers metadata if QUARC used (engine should set this flag in result later in Phase 1)

    # If --output-file is provided, write UTF-8 to that file as well
    out_path = _extract_output_path()
    if out_path:
        try:
            with open(out_path, 'w', encoding='utf-8') as f:
                json.dump(export, f, ensure_ascii=False)
        except Exception:
            # ignore file write errors; still write to stdout
            pass

    _write_json(export)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
