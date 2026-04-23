#!/usr/bin/env python3
"""CLIP-seq skill — mines flow.bio execution history to suggest parameters for the next CLIP-seq run."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

CLIP_PIPELINE_ID = "960154035051242353"
CLIP_PIPELINE_VERSION_ID = "337265464908053502"  # v1.7 (latest)
CLIP_PIPELINE_NAME = "CLIP-Seq"
VERSION = "0.2.0"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

# Schema for all tunable pipeline parameters (from flow.bio pipeline version schema)
PARAM_SCHEMA: dict[str, dict] = {
    "crosslink_position": {
        "type": "categorical",
        "choices": ["start", "end", "middle"],
        "default": "start",
        "description": "Position of the crosslink in the read. 'start' for iCLIP, 'end' for eCLIP.",
    },
    "encode_eclip": {
        "type": "boolean",
        "default": False,
        "description": "Enable ENCODE eCLIP adapter/UMI settings.",
    },
    "move_umi_to_header": {
        "type": "boolean",
        "default": False,
        "description": "Move inline UMI to read header before alignment.",
    },
    "umi_separator": {
        "type": "string",
        "default": "rbc:",
        "description": "Delimiter UMICollapse uses to find UMIs in aligned read names.",
    },
    "skip_umi_dedupe": {
        "type": "boolean",
        "default": False,
        "description": "Skip UMI-based PCR duplicate removal.",
    },
    "paraclu_min_value": {
        "type": "number",
        "default": None,
        "description": "Minimum cluster value for Paraclu peak calling.",
    },
    "trimgalore_params": {
        "type": "string",
        "default": None,
        "description": "Extra arguments passed to TrimGalore.",
    },
    "star_params": {
        "type": "string",
        "default": None,
        "description": "Extra arguments passed to STAR alignment.",
    },
    "bowtie_params": {
        "type": "string",
        "default": None,
        "description": "Extra arguments passed to Bowtie pre-mapping (smRNA removal).",
    },
    "clippy_params": {
        "type": "string",
        "default": None,
        "description": "Extra arguments passed to Clippy peak caller.",
    },
    "icount_peak_params": {
        "type": "string",
        "default": None,
        "description": "Extra arguments passed to iCount peak calling.",
    },
    "peka_params": {
        "type": "string",
        "default": None,
        "description": "Extra arguments passed to PEKA K-mer enrichment analysis.",
    },
}

# ---------------------------------------------------------------------------
# Core analysis functions
# ---------------------------------------------------------------------------


def extract_params(execution: dict) -> dict:
    """Extract and normalise tunable parameters from a single execution dict.

    Coerces string booleans ('true'/'false') to Python bools.
    Adds 'genome' from the fileset field.
    """
    raw = execution.get("params") or {}
    result: dict = {}

    for key, value in raw.items():
        if isinstance(value, str):
            if value.lower() == "true":
                value = True
            elif value.lower() == "false":
                value = False
        result[key] = value

    fileset = execution.get("fileset") or {}
    result["genome"] = fileset.get("name") or fileset.get("organism_name") or "unknown" if fileset else "unknown"

    return result


def aggregate_params(param_list: list[dict]) -> dict[str, dict]:
    """Count occurrences of each value for every parameter across all runs.

    Returns a dict mapping param_name → {value: count, ...}.
    Only includes parameters that appear in at least one execution.
    """
    counters: dict[str, Counter] = {}
    for params in param_list:
        for key, value in params.items():
            if key not in counters:
                counters[key] = Counter()
            counters[key][value] += 1

    return {k: dict(v) for k, v in counters.items()}


def suggest_params(stats: dict[str, dict]) -> dict[str, dict]:
    """Derive a suggested value + confidence for each parameter.

    Confidence rules:
      - 1 run total  → 'low'
      - 2–4 runs     → 'medium'
      - 5+ runs with ≥ 70% agreement → 'high'
      - 5+ runs with < 70% agreement → 'medium'
    """
    suggestions: dict[str, dict] = {}

    for param, value_counts in stats.items():
        total = sum(value_counts.values())
        most_common_value = max(value_counts, key=value_counts.__getitem__)
        top_count = value_counts[most_common_value]
        agreement = top_count / total if total > 0 else 0

        if total == 1:
            confidence = "low"
        elif total < 5:
            confidence = "medium"
        elif agreement >= 0.70:
            confidence = "high"
        else:
            confidence = "medium"

        suggestions[param] = {
            "value": most_common_value,
            "confidence": confidence,
            "runs_analysed": total,
            "agreement_pct": round(agreement * 100, 1),
            "frequency_table": value_counts,
        }

    return suggestions


# ---------------------------------------------------------------------------
# Demo data
# ---------------------------------------------------------------------------


def _make_demo_executions() -> list[dict]:
    """Return a list of synthetic CLIP-seq execution dicts for demo/testing.

    Represents a realistic mix of iCLIP and eCLIP parameter choices across
    15 runs from 3 different labs.
    """
    base_ts = 1_700_000_000

    runs = [
        # iCLIP runs (typical: crosslink_position=start, no ENCODE eCLIP)
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",
         "star_params": "--outFilterMismatchNmax 2", "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "true",
         "move_umi_to_header": "false", "umi_separator": "rbc:", "status": "ERROR"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",
         "trimgalore_params": "--clip_R1 5", "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
        # eCLIP runs (crosslink_position=end, encode_eclip=true)
        {"crosslink_position": "end",   "encode_eclip": "true",  "skip_umi_dedupe": "false",
         "move_umi_to_header": "false", "umi_separator": "rbc:", "status": "OK"},
        {"crosslink_position": "end",   "encode_eclip": "true",  "skip_umi_dedupe": "false",
         "move_umi_to_header": "false", "umi_separator": "rbc:", "status": "OK"},
        {"crosslink_position": "end",   "encode_eclip": "true",  "skip_umi_dedupe": "false",
         "move_umi_to_header": "false", "umi_separator": "rbc:",
         "star_params": "--outFilterMismatchNmax 2", "status": "OK"},
        # Mixed / exploratory runs
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",
         "paraclu_min_value": "30",     "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",
         "paraclu_min_value": "50",     "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
        {"crosslink_position": "start", "encode_eclip": "false", "skip_umi_dedupe": "false",
         "move_umi_to_header": "true",  "umi_separator": "_",    "status": "OK"},
    ]

    executions = []
    for i, run in enumerate(runs):
        params = {k: v for k, v in run.items() if k != "status"}
        executions.append({
            "id": str(800_000_000_000_000_000 + i),
            "identifier": f"demo_run_{i+1:02d}",
            "pipeline_name": CLIP_PIPELINE_NAME,
            "pipeline_version": "1.7",
            "status": run["status"],
            "created": base_ts + i * 86400,
            "params": params,
            "fileset": {"name": "GRCh38", "organism_name": "Human"},
        })

    return executions


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _format_frequency_table(freq: dict) -> str:
    """Format a value-count dict as a markdown table row string."""
    parts = [f"`{v}` × {c}" for v, c in sorted(freq.items(), key=lambda x: -x[1])]
    return ", ".join(parts)


def generate_report(
    executions: list[dict],
    suggestions: dict[str, dict],
    stats: dict[str, dict],
    output_dir: Path,
    mode: str = "demo",
) -> Path:
    n = len(executions)
    ok = sum(1 for e in executions if e.get("status") == "OK")
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = [
        "# ClawBio CLIP-seq Parameter Advisor",
        "",
        f"**Date**: {ts}",
        f"**Mode**: {mode}",
        f"**Executions analysed**: {n} ({ok} successful)",
        f"**Pipeline**: {CLIP_PIPELINE_NAME} on flow.bio",
        "",
        "---",
        "",
        "## Suggested Parameters for Next Run",
        "",
        "Based on the execution history above, the following parameters are recommended:",
        "",
        "| Parameter | Suggested Value | Confidence | Runs | Agreement |",
        "|-----------|----------------|------------|------|-----------|",
    ]

    for param, s in suggestions.items():
        val = s["value"]
        conf_emoji = {"high": "✅", "medium": "⚠️", "low": "❓"}.get(s["confidence"], "")
        lines.append(
            f"| `{param}` | `{val}` | {conf_emoji} {s['confidence']} "
            f"| {s['runs_analysed']} | {s['agreement_pct']}% |"
        )

    lines += [
        "",
        "## Parameter Frequency Distribution",
        "",
    ]

    for param, freq in stats.items():
        schema_info = PARAM_SCHEMA.get(param, {})
        desc = schema_info.get("description", "")
        lines.append(f"### `{param}`")
        if desc:
            lines.append(f"*{desc}*")
        lines.append("")
        lines.append("| Value | Count |")
        lines.append("|-------|-------|")
        for val, count in sorted(freq.items(), key=lambda x: -x[1]):
            lines.append(f"| `{val}` | {count} |")
        lines.append("")

    lines += [
        "## Execution History",
        "",
        "| # | ID | Status | crosslink_position | encode_eclip | skip_umi_dedupe | genome |",
        "|---|-----|--------|--------------------|--------------|-----------------|--------|",
    ]
    for i, e in enumerate(executions, 1):
        p = e.get("params", {})
        lines.append(
            f"| {i} | `{str(e.get('id','?'))[:12]}…` | {e.get('status','?')} "
            f"| {p.get('crosslink_position','—')} "
            f"| {p.get('encode_eclip','—')} "
            f"| {p.get('skip_umi_dedupe','—')} "
            f"| {(e.get('fileset') or {}).get('name','—')} |"
        )

    lines += [
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]

    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


# ---------------------------------------------------------------------------
# Run modes
# ---------------------------------------------------------------------------


def run_demo(output_dir: Path) -> None:
    """Run with synthetic execution history — no credentials required."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    executions = _make_demo_executions()
    param_list = [extract_params(e) for e in executions]
    stats = aggregate_params(param_list)
    suggestions = suggest_params(stats)

    generate_report(executions, suggestions, stats, output_dir, mode="DEMO")

    result = {
        "skill": "clip-seq",
        "version": VERSION,
        "mode": "demo",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "executions_analysed": executions,
        "param_stats": stats,
        "suggestions": {k: {"value": v["value"], "confidence": v["confidence"]} for k, v in suggestions.items()},
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8"
    )

    (output_dir / "reproducibility" / "commands.sh").write_text(
        f"#!/usr/bin/env bash\n"
        f"# ClawBio CLIP-seq demo — {datetime.now(timezone.utc).isoformat()}\n\n"
        f"python skills/clip-seq/clip_seq.py --demo --output /tmp/clip_seq_demo\n",
        encoding="utf-8",
    )

    print(f"Demo report written to {output_dir}/report.md")
    print(f"Suggestions based on {len(executions)} synthetic executions.")


def run_history(
    username: str,
    password: str,
    output_dir: Path,
    n_executions: int = 50,
    include_public: bool = True,
) -> None:
    """Fetch real CLIP-seq execution history from flow.bio and analyse parameters."""
    try:
        sys.path.insert(0, str(SKILL_DIR.parent / "flow-bio"))
        from flow_bio import FlowClient
    except ImportError:
        print("ERROR: flow-bio skill not found. Run from the ClawBio root directory.", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    client = FlowClient()
    resp = client.login(username, password)
    user = resp.get("user", {}).get("username", username)
    print(f"Authenticated as: {user}")

    # Fetch owned executions and filter to CLIP-Seq
    print("Fetching CLIP-Seq execution history...")
    owned = client.get_executions_owned()
    clip_owned = [e for e in owned if CLIP_PIPELINE_NAME in e.get("pipeline_name", "")]
    print(f"  Owned CLIP-Seq executions: {len(clip_owned)}")

    # If too few owned runs, supplement with public executions via search
    executions = list(clip_owned)
    if include_public and len(executions) < n_executions:
        needed = n_executions - len(executions)
        print(f"  Fetching up to {needed} public CLIP-Seq executions for richer parameter analysis...")
        try:
            # Use a known public CLIP sample to find executions
            pub_result = client._get("/executions/search", params={"pipeline_name": CLIP_PIPELINE_NAME})
            pub_execs = pub_result.get("executions", []) if isinstance(pub_result, dict) else []
            # Fetch detail for a sample of public runs to get params
            for e in pub_execs[:needed]:
                try:
                    detail = client.get_execution(e["id"])
                    if detail.get("params"):
                        executions.append(detail)
                except Exception:
                    pass
        except Exception as exc:
            print(f"  Warning: could not fetch public executions ({exc})", file=sys.stderr)

    if not executions:
        print("No CLIP-Seq executions found. Run some pipelines first, or use --demo.", file=sys.stderr)
        sys.exit(1)

    # For owned executions that only have summary fields, fetch details
    detailed: list[dict] = []
    for e in executions:
        if e.get("params") is None:
            try:
                detailed.append(client.get_execution(e["id"]))
            except Exception:
                detailed.append(e)
        else:
            detailed.append(e)

    param_list = [extract_params(e) for e in detailed]
    stats = aggregate_params(param_list)
    suggestions = suggest_params(stats)

    generate_report(detailed, suggestions, stats, output_dir, mode="live")

    result = {
        "skill": "clip-seq",
        "version": VERSION,
        "mode": "live",
        "user": user,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "executions_analysed": detailed,
        "param_stats": {k: {str(vv): c for vv, c in v.items()} for k, v in stats.items()},
        "suggestions": {k: {"value": v["value"], "confidence": v["confidence"]} for k, v in suggestions.items()},
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8"
    )

    (output_dir / "reproducibility" / "commands.sh").write_text(
        f"#!/usr/bin/env bash\n"
        f"# ClawBio CLIP-seq history — {datetime.now(timezone.utc).isoformat()}\n\n"
        f"FLOW_USERNAME='{username}' FLOW_PASSWORD='***' \\\n"
        f"  python skills/clip-seq/clip_seq.py --history --output {output_dir}\n",
        encoding="utf-8",
    )

    print(f"\nReport: {output_dir}/report.md")
    print(f"Analysed {len(detailed)} executions.")
    print("\nTop suggestions:")
    for param, s in list(suggestions.items())[:5]:
        print(f"  {param}: {s['value']}  ({s['confidence']} confidence, {s['agreement_pct']}% agreement)")


# ---------------------------------------------------------------------------
# --run mode helpers
# ---------------------------------------------------------------------------


def build_pipeline_params(
    suggestions: dict[str, dict],
    override: dict | None = None,
) -> dict[str, str]:
    """Convert suggestion dict to flat string params ready for flow.bio API.

    - Skips 'genome' (passed as fileset ID separately).
    - Skips params not in PARAM_SCHEMA (e.g. internal keys).
    - Skips None values (unset optional params).
    - Coerces booleans to 'true'/'false' strings.
    - Override values take precedence over suggestions.
    """
    params: dict[str, str] = {}

    for param, s in suggestions.items():
        if param == "genome" or param not in PARAM_SCHEMA:
            continue
        val = s["value"]
        if val is None:
            continue
        params[param] = "true" if val is True else ("false" if val is False else str(val))

    if override:
        for k, v in override.items():
            if v is None:
                continue
            params[k] = "true" if v is True else ("false" if v is False else str(v))

    return params


def run_pipeline_with_suggestions(
    username: str,
    password: str,
    sample_id: str,
    genome_id: str,
    output_dir: Path,
    pipeline_version_id: str = CLIP_PIPELINE_VERSION_ID,
    n_history: int = 50,
    override_params: dict | None = None,
    dry_run: bool = False,
) -> None:
    """Mine execution history, build params, and launch (or preview) the pipeline.

    In dry-run mode: shows what would be submitted, writes report + result.json,
    exits without calling the flow.bio run endpoint.
    Falls back to demo suggestions if authentication fails (useful for dry-run testing).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    # --- Fetch history for suggestions ---
    executions: list[dict] = []
    suggestions: dict[str, dict] = {}
    history_source = "demo fallback"

    try:
        sys.path.insert(0, str(SKILL_DIR.parent / "flow-bio"))
        from flow_bio import FlowClient

        client = FlowClient()
        resp = client.login(username, password)
        user = resp.get("user", {}).get("username", username)
        print(f"Authenticated as: {user}")

        owned = client.get_executions_owned()
        clip_owned = [e for e in owned if CLIP_PIPELINE_NAME in e.get("pipeline_name", "")]

        executions = list(clip_owned)
        if len(executions) < n_history:
            try:
                pub = client._get("/executions/search", params={"pipeline_name": CLIP_PIPELINE_NAME})
                for e in (pub.get("executions", []) if isinstance(pub, dict) else [])[:n_history - len(executions)]:
                    try:
                        executions.append(client.get_execution(e["id"]))
                    except Exception:
                        pass
            except Exception:
                pass

        # Fetch params for owned summary-only records
        detailed: list[dict] = []
        for e in executions:
            if e.get("params") is None:
                try:
                    detailed.append(client.get_execution(e["id"]))
                except Exception:
                    detailed.append(e)
            else:
                detailed.append(e)

        executions = detailed
        history_source = f"flow.bio ({len(executions)} executions)"

    except Exception as exc:
        if not dry_run:
            print(f"ERROR: authentication failed — {exc}", file=sys.stderr)
            sys.exit(1)
        print(f"Warning: could not authenticate ({exc}). Using demo suggestions for dry run.", file=sys.stderr)
        executions = _make_demo_executions()
        history_source = "demo data (auth unavailable)"

    param_list = [extract_params(e) for e in executions]
    stats = aggregate_params(param_list)
    suggestions = suggest_params(stats)
    pipeline_params = build_pipeline_params(suggestions, override=override_params)

    # --- Report ---
    mode_label = f"{'DRY RUN — ' if dry_run else ''}params from {history_source}"
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = [
        "# ClawBio CLIP-seq — Pipeline Launch" + (" (DRY RUN)" if dry_run else ""),
        "",
        f"**Date**: {ts}",
        f"**Sample ID**: `{sample_id}`",
        f"**Genome fileset ID**: `{genome_id}`",
        f"**Pipeline version**: `{pipeline_version_id}`",
        f"**Mode**: {mode_label}",
        "",
    ]

    if dry_run:
        lines += [
            "> **DRY RUN** — No pipeline was launched. "
            "Remove `--dry-run` to submit.",
            "",
        ]

    lines += [
        "## Parameters to be submitted",
        "",
        "| Parameter | Value | Source |",
        "|-----------|-------|--------|",
    ]
    for k, v in pipeline_params.items():
        src = "override" if (override_params and k in override_params) else \
              f"{suggestions.get(k, {}).get('confidence', '?')} confidence"
        lines.append(f"| `{k}` | `{v}` | {src} |")

    lines += [
        "",
        "## Suggested values (full detail)",
        "",
        "| Parameter | Value | Confidence | Agreement |",
        "|-----------|-------|------------|-----------|",
    ]
    for param, s in suggestions.items():
        if param == "genome" or param not in PARAM_SCHEMA:
            continue
        conf_emoji = {"high": "✅", "medium": "⚠️", "low": "❓"}.get(s["confidence"], "")
        lines.append(
            f"| `{param}` | `{s['value']}` | {conf_emoji} {s['confidence']} | {s['agreement_pct']}% |"
        )

    if not dry_run:
        lines += ["", f"## Execution", "", f"**Execution ID**: see `result.json`"]

    lines += [
        "",
        "## Reproduce",
        "",
        "```bash",
        f"python skills/clip-seq/clip_seq.py \\",
        f"  --run --sample {sample_id} --genome {genome_id} \\",
        f"  --output {output_dir}",
        "```",
        "",
        f"---",
        f"*{DISCLAIMER}*",
        "",
    ]

    (output_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # --- Launch or dry-run ---
    execution_result: dict = {}
    if not dry_run:
        try:
            from flow_bio import FlowClient as _FC
            _client = _FC()
            _client.login(username, password)
            execution_result = _client.run_pipeline(
                pipeline_version_id=pipeline_version_id,
                sample_ids=[sample_id],
                params=pipeline_params,
                genome_id=genome_id,
            )
            exec_id = execution_result.get("id", "?")
            print(f"\nPipeline launched! Execution ID: {exec_id}")
            print(f"Monitor with: python skills/flow-bio/flow_bio.py --execution {exec_id}")
        except Exception as exc:
            print(f"ERROR: pipeline launch failed — {exc}", file=sys.stderr)
            sys.exit(1)

    result = {
        "skill": "clip-seq",
        "version": VERSION,
        "mode": "dry_run" if dry_run else "run",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sample_id": sample_id,
        "genome_id": genome_id,
        "pipeline_version_id": pipeline_version_id,
        "pipeline_params": pipeline_params,
        "suggestions": {k: {"value": v["value"], "confidence": v["confidence"]} for k, v in suggestions.items()},
        "execution": execution_result,
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8"
    )

    cmd = (
        f"FLOW_USERNAME='{username}' FLOW_PASSWORD='***' \\\n"
        f"  python skills/clip-seq/clip_seq.py \\\n"
        f"  --run --sample {sample_id} --genome {genome_id} \\\n"
        f"  --output {output_dir}\n"
    )
    (output_dir / "reproducibility" / "commands.sh").write_text(
        f"#!/usr/bin/env bash\n# {ts}\n\n{cmd}", encoding="utf-8"
    )

    if dry_run:
        print(f"\nDry run complete. Review params in {output_dir}/report.md before launching.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="CLIP-seq: mine flow.bio execution history to suggest pipeline parameters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python clip_seq.py --demo --output /tmp/clip_demo\n"
            "  python clip_seq.py --history --output /tmp/clip_history\n"
            "  python clip_seq.py --run --sample <ID> --genome <ID> --dry-run --output /tmp/out\n"
            "  python clip_seq.py --run --sample <ID> --genome <ID> --output /tmp/out\n"
        ),
    )
    p.add_argument("--demo",    action="store_true", help="Run with synthetic data (no credentials)")
    p.add_argument("--history", action="store_true", help="Fetch CLIP-Seq history and show suggestions")
    p.add_argument("--run",     action="store_true", help="Fetch history, build params, and launch pipeline")

    p.add_argument("--username", help="flow.bio username (or set FLOW_USERNAME)")
    p.add_argument("--password", help="flow.bio password (or set FLOW_PASSWORD)")
    p.add_argument("--n",        type=int, default=50, help="Max executions to analyse (default: 50)")
    p.add_argument("--no-public", action="store_true", help="Only use owned executions")

    # --run specific
    p.add_argument("--sample",   metavar="ID",   help="flow.bio sample ID to run the pipeline on")
    p.add_argument("--genome",   metavar="ID",   help="flow.bio genome fileset ID")
    p.add_argument("--override", metavar="JSON", help="JSON string of params that override suggestions")
    p.add_argument("--dry-run",  action="store_true", help="Show params without launching the pipeline")
    p.add_argument("--pipeline-version", default=CLIP_PIPELINE_VERSION_ID,
                   help=f"Pipeline version ID (default: {CLIP_PIPELINE_VERSION_ID})")

    p.add_argument("--output", required=True, help="Output directory")
    return p.parse_args()


def main() -> None:
    import os
    args = parse_args()
    output_dir = Path(args.output)

    if args.demo:
        run_demo(output_dir)
        return

    username = args.username or os.environ.get("FLOW_USERNAME", "")
    password = args.password or os.environ.get("FLOW_PASSWORD", "")

    if args.history:
        if not username or not password:
            print(
                "ERROR: --history requires flow.bio credentials.\n"
                "Set FLOW_USERNAME + FLOW_PASSWORD, or use --username / --password.\n"
                "Try --demo to run without credentials.",
                file=sys.stderr,
            )
            sys.exit(1)
        run_history(
            username=username,
            password=password,
            output_dir=output_dir,
            n_executions=args.n,
            include_public=not args.no_public,
        )
        return

    if args.run:
        if not args.sample or not args.genome:
            print("ERROR: --run requires --sample <ID> and --genome <ID>", file=sys.stderr)
            sys.exit(1)
        override = json.loads(args.override) if args.override else None
        run_pipeline_with_suggestions(
            username=username,
            password=password,
            sample_id=args.sample,
            genome_id=args.genome,
            output_dir=output_dir,
            pipeline_version_id=args.pipeline_version,
            n_history=args.n,
            override_params=override,
            dry_run=args.dry_run,
        )
        return

    print("ERROR: specify --demo, --history, or --run", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
