"""Tests for the clip-seq skill (red/green TDD — fill these in before implementing)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).parent.parent
SCRIPT = SKILL_DIR / "clip_seq.py"


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

def test_demo_runs(tmp_path):
    """Demo mode exits 0 and writes report.md."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "report.md").exists()


def test_demo_report_contains_disclaimer(tmp_path):
    """Demo report includes the ClawBio medical disclaimer."""
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    text = (tmp_path / "report.md").read_text()
    assert "research and educational tool" in text.lower() or "not a medical device" in text.lower()


def test_demo_creates_output_structure(tmp_path):
    """Demo creates expected subdirectories."""
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    assert (tmp_path / "figures").is_dir()
    assert (tmp_path / "tables").is_dir()
    assert (tmp_path / "reproducibility").is_dir()


def test_demo_metadata_json(tmp_path):
    """Demo writes valid metadata JSON with skill name and version."""
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True,
    )
    meta_path = tmp_path / "tables" / "metadata.json"
    assert meta_path.exists()
    meta = json.loads(meta_path.read_text())
    assert meta["skill"] == "clip-seq"
    assert "version" in meta


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def test_missing_input_exits_nonzero(tmp_path):
    """Running without --input or --demo should exit with an error."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0


def test_nonexistent_input_exits_nonzero(tmp_path):
    """Passing a file that does not exist should exit with an error."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", "/nonexistent/file.fastq.gz",
         "--output", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0


# ---------------------------------------------------------------------------
# Full pipeline (add tests here as you implement each step)
# ---------------------------------------------------------------------------

# TODO: test adapter trimming output
# TODO: test alignment produces sorted BAM
# TODO: test peak calling produces non-empty BED
# TODO: test peak annotation maps to expected features
# TODO: test UMI deduplication reduces read count
