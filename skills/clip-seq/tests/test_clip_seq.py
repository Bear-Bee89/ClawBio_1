"""Tests for the clip-seq skill (red/green TDD)."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).parent.parent
SCRIPT = SKILL_DIR / "clip_seq.py"

CLIP_PIPELINE_ID = "960154035051242353"


def _load_module():
    spec = importlib.util.spec_from_file_location("clip_seq", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Demo mode — CLI
# ---------------------------------------------------------------------------

def test_demo_runs(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr


def test_demo_creates_output_structure(tmp_path):
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "result.json").exists()
    assert (tmp_path / "reproducibility").is_dir()


def test_demo_report_contains_disclaimer(tmp_path):
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    text = (tmp_path / "report.md").read_text()
    assert "not a medical device" in text.lower() or "research and educational tool" in text.lower()


def test_demo_report_has_param_suggestions(tmp_path):
    """Report must contain a parameter suggestions section."""
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    text = (tmp_path / "report.md").read_text()
    assert "crosslink_position" in text


def test_demo_result_json_has_suggestions(tmp_path):
    """result.json must contain a 'suggestions' key with at least one param."""
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    data = json.loads((tmp_path / "result.json").read_text())
    assert "suggestions" in data
    assert "crosslink_position" in data["suggestions"]


def test_demo_result_json_has_execution_history(tmp_path):
    """result.json must include the executions list used for analysis."""
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    data = json.loads((tmp_path / "result.json").read_text())
    assert "executions_analysed" in data
    assert len(data["executions_analysed"]) >= 5


# ---------------------------------------------------------------------------
# Unit tests — importable functions
# ---------------------------------------------------------------------------

def test_extract_params_returns_dict():
    """extract_params must return a flat dict of pipeline parameter values."""
    mod = _load_module()
    execution = {
        "id": "abc123",
        "status": "OK",
        "pipeline_name": "CLIP-Seq",
        "params": {
            "crosslink_position": "start",
            "encode_eclip": "false",
            "skip_umi_dedupe": "false",
            "move_umi_to_header": "false",
            "umi_separator": "_",
        },
        "fileset": {"name": "GRCh38", "organism_name": "Human"},
    }
    params = mod.extract_params(execution)
    assert isinstance(params, dict)
    assert params["crosslink_position"] == "start"
    assert params["genome"] == "GRCh38"


def test_extract_params_coerces_booleans():
    """String 'true'/'false' should be coerced to Python bools."""
    mod = _load_module()
    execution = {
        "params": {"skip_umi_dedupe": "true", "encode_eclip": "false"},
        "fileset": {},
    }
    params = mod.extract_params(execution)
    assert params["skip_umi_dedupe"] is True
    assert params["encode_eclip"] is False


def test_extract_params_handles_missing_fileset():
    """extract_params must not raise if fileset is absent."""
    mod = _load_module()
    execution = {
        "params": {"crosslink_position": "end"},
    }
    params = mod.extract_params(execution)
    assert params.get("genome") in (None, "unknown")


def test_aggregate_params_counts_values():
    """aggregate_params must count occurrences of each value per parameter."""
    mod = _load_module()
    param_list = [
        {"crosslink_position": "start", "skip_umi_dedupe": False},
        {"crosslink_position": "start", "skip_umi_dedupe": True},
        {"crosslink_position": "end",   "skip_umi_dedupe": False},
    ]
    stats = mod.aggregate_params(param_list)
    assert stats["crosslink_position"]["start"] == 2
    assert stats["crosslink_position"]["end"] == 1
    assert stats["skip_umi_dedupe"][False] == 2
    assert stats["skip_umi_dedupe"][True] == 1


def test_suggest_params_picks_most_common():
    """suggest_params must recommend the most frequently used value."""
    mod = _load_module()
    stats = {
        "crosslink_position": {"start": 8, "end": 2},
        "encode_eclip": {False: 9, True: 1},
    }
    suggestions = mod.suggest_params(stats)
    assert suggestions["crosslink_position"]["value"] == "start"
    assert suggestions["encode_eclip"]["value"] is False


def test_suggest_params_includes_confidence():
    """Each suggestion must include a confidence level."""
    mod = _load_module()
    stats = {"crosslink_position": {"start": 10}}
    suggestions = mod.suggest_params(stats)
    assert "confidence" in suggestions["crosslink_position"]
    assert suggestions["crosslink_position"]["confidence"] in ("high", "medium", "low")


def test_suggest_params_low_confidence_for_single_run():
    """A param seen in only one run gets 'low' confidence."""
    mod = _load_module()
    stats = {"crosslink_position": {"end": 1}}
    suggestions = mod.suggest_params(stats)
    assert suggestions["crosslink_position"]["confidence"] == "low"


def test_make_demo_executions_returns_list():
    """_make_demo_executions must return a non-empty list of execution dicts."""
    mod = _load_module()
    execs = mod._make_demo_executions()
    assert isinstance(execs, list)
    assert len(execs) >= 5
    for e in execs:
        assert "params" in e
        assert "status" in e


def test_make_demo_executions_has_varied_params():
    """Demo executions should not all have identical parameters."""
    mod = _load_module()
    execs = mod._make_demo_executions()
    positions = {e["params"].get("crosslink_position") for e in execs if e.get("params")}
    assert len(positions) >= 1


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def test_missing_output_exits_nonzero():
    """Running with no --output or --demo should exit non-zero."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0


# ---------------------------------------------------------------------------
# Full pipeline — add tests as each step is implemented
# ---------------------------------------------------------------------------

# TODO: test --history mode authenticates and fetches real executions
# TODO: test param extraction from a live execution
# TODO: test report reflects actual execution count
