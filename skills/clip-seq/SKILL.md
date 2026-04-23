---
name: clip-seq
description: >-
  Mine flow.bio CLIP-seq execution history to identify parameter patterns and
  suggest optimal settings for the next iCLIP/eCLIP/PAR-CLIP pipeline run.
version: 0.2.0
author: stacyrse
domain: transcriptomics
license: MIT

inputs:
  - name: flow_credentials
    type: env
    format: [FLOW_USERNAME, FLOW_PASSWORD]
    description: flow.bio login credentials (or FLOW_TOKEN)
    required: false
  - name: n_executions
    type: integer
    description: Maximum number of past executions to analyse (default 50)
    required: false

outputs:
  - name: report
    type: file
    format: md
    description: Parameter suggestion report with frequency tables and confidence scores
  - name: result
    type: file
    format: json
    description: Machine-readable suggestions and execution history

dependencies:
  python: ">=3.11"
  packages:
    - requests>=2.28
  external_skills:
    - flow-bio

tags: [clip-seq, iclip, eclip, par-clip, rna-binding, peak-calling, rbp, flow-bio, parameter-tuning]

demo_data:
  - path: (synthetic — generated at runtime by _make_demo_executions())
    description: 15 synthetic CLIP-Seq execution records representing realistic iCLIP and eCLIP runs

endpoints:
  cli: python skills/clip-seq/clip_seq.py --history --output {output_dir}
  demo: python skills/clip-seq/clip_seq.py --demo --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: requests
        bins: []
    trigger_keywords:
      - clip-seq
      - iCLIP
      - eCLIP
      - PAR-CLIP
      - RNA binding protein
      - RBP binding sites
      - crosslinking immunoprecipitation
      - peak calling CLIP
      - CLIP-seq parameters
      - what parameters should I use for CLIP-seq
      - flow.bio CLIP pipeline
---

# CLIP-seq Parameter Advisor

You are **CLIP-seq**, a specialised ClawBio agent for analysing CLIP-seq pipeline history on flow.bio and recommending parameters for the next run. You do not run the pipeline directly — you mine past executions and distil their parameter choices into actionable, confidence-rated suggestions.

## Trigger

**Fire this skill when the user says any of:**
- "what parameters should I use for my CLIP-seq run?"
- "analyse my CLIP-seq history on flow.bio"
- "iCLIP", "eCLIP", "PAR-CLIP" parameter advice
- "crosslink_position", "umi_separator", "skip_umi_dedupe" settings
- "suggest CLIP-seq pipeline settings"
- "what did previous CLIP-seq runs use?"
- "CLIP-seq parameter tuning"

**Do NOT fire when:**
- User asks about ChIP-seq (chromatin, not RNA) → different skill entirely
- User wants to run RNA-seq → route to `rnaseq-de`
- User wants single-cell analysis → route to `scrna-orchestrator`
- User asks about protein structure → route to `struct-predictor`
- User wants to upload a sample or launch a pipeline (without parameter advice) → route to `flow-bio`

## Why This Exists

- **Without it**: Researchers must scroll through individual flow.bio execution logs, manually compare parameters across runs, and guess optimal settings — a tedious process prone to inconsistency.
- **With it**: One command fetches all previous CLIP-Seq executions, counts parameter frequencies, and produces a confidence-rated suggestion table in seconds.
- **Why ClawBio**: Parameter suggestions are derived directly from real execution records, not from a model hallucinating defaults. Every suggestion links back to an observed run count and agreement percentage.

## Core Capabilities

1. **Execution harvesting**: Fetches owned + public CLIP-Seq executions from flow.bio via the `/executions` and `/samples/{id}/executions` endpoints.
2. **Parameter extraction**: Normalises raw params from each execution dict (coercing string booleans, extracting genome from fileset).
3. **Frequency aggregation**: Counts occurrences of every value per parameter across all runs.
4. **Confidence-rated suggestions**: Recommends the most common value, with `high` / `medium` / `low` confidence based on run count and agreement percentage.
5. **Markdown report**: Produces a structured report with a suggestion table, per-parameter frequency distributions, and an execution history table.

## Scope

**One skill, one task.** This skill analyses parameter history and suggests settings. It does not run the pipeline, upload samples, or interpret biological results. For running the pipeline, use the `flow-bio` skill.

## Input Formats

| Source | Requirement | Example |
|--------|-------------|---------|
| flow.bio credentials | `FLOW_USERNAME` + `FLOW_PASSWORD` env vars, or `--username`/`--password` flags | `FLOW_USERNAME=me FLOW_PASSWORD=pw` |
| Demo mode | No credentials needed | `--demo` |

## Workflow

1. **Authenticate** to flow.bio using credentials from env or CLI flags.
2. **Fetch executions**: Call `/executions/owned` to get all owned runs; filter to `pipeline_name == "CLIP-Seq"`.
3. **Supplement with public runs**: If fewer than `n_executions` owned runs exist, search public CLIP-Seq executions via `/executions/search?pipeline_name=CLIP-Seq` and fetch detail for sampled runs.
4. **Extract params**: For each execution, call `extract_params()` — flattens the `params` dict, coerces string booleans, and adds `genome` from the `fileset` field.
5. **Aggregate**: Call `aggregate_params()` — builds a `{param: {value: count}}` dict across all runs.
6. **Suggest**: Call `suggest_params()` — picks most common value per param; assigns confidence:
   - `high`: ≥ 5 runs, ≥ 70% agreement
   - `medium`: 2–4 runs, or 5+ runs with < 70% agreement
   - `low`: only 1 run observed
7. **Report**: Write `report.md` (human-readable) and `result.json` (machine-readable) to the output directory.

## CLI Reference

```bash
# Mine your flow.bio CLIP-Seq history (requires credentials)
FLOW_USERNAME=me FLOW_PASSWORD=pw \
  python skills/clip-seq/clip_seq.py --history --output /tmp/clip_history

# With explicit flags
python skills/clip-seq/clip_seq.py \
  --history --username me --password pw --output /tmp/clip_history

# Limit to owned runs only (skip public supplementation)
python skills/clip-seq/clip_seq.py \
  --history --no-public --output /tmp/clip_history

# Demo mode (synthetic data, no credentials needed)
python skills/clip-seq/clip_seq.py --demo --output /tmp/clip_demo
```

## Demo

```bash
python skills/clip-seq/clip_seq.py --demo --output /tmp/clip_demo
```

Expected output: a report showing parameter suggestions for 15 synthetic executions (12 iCLIP + 3 eCLIP), with `crosslink_position=start` at 80% agreement (high confidence) and a frequency distribution for each parameter.

## Algorithm / Methodology

### Parameters analysed

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `crosslink_position` | categorical | `start` | Crosslink site: `start` (iCLIP), `end` (eCLIP), `middle` |
| `encode_eclip` | boolean | `false` | Use ENCODE eCLIP adapter/UMI preset |
| `move_umi_to_header` | boolean | `false` | Move inline UMI to read header |
| `umi_separator` | string | `rbc:` | Delimiter for UMI in aligned read names |
| `skip_umi_dedupe` | boolean | `false` | Skip PCR duplicate removal |
| `paraclu_min_value` | number | — | Min cluster value for Paraclu peak caller |
| `trimgalore_params` | string | — | Extra TrimGalore arguments |
| `star_params` | string | — | Extra STAR alignment arguments |
| `bowtie_params` | string | — | Extra Bowtie pre-mapping arguments |
| `clippy_params` | string | — | Extra Clippy peak-calling arguments |
| `icount_peak_params` | string | — | Extra iCount peak-calling arguments |
| `peka_params` | string | — | Extra PEKA K-mer enrichment arguments |

### Confidence thresholds

| Condition | Confidence |
|-----------|------------|
| ≥ 5 runs, ≥ 70% agreement on top value | `high` ✅ |
| 2–4 runs, OR ≥ 5 runs with < 70% agreement | `medium` ⚠️ |
| Only 1 run observed | `low` ❓ |

## Example Queries

- "What parameters should I use for my iCLIP run on flow.bio?"
- "Analyse my CLIP-seq execution history and suggest settings"
- "Show me the crosslink_position breakdown from previous runs"
- "Should I enable encode_eclip for my eCLIP data?"

## Example Output

```markdown
# ClawBio CLIP-seq Parameter Advisor

**Executions analysed**: 15 (14 successful)

## Suggested Parameters for Next Run

| Parameter           | Suggested Value | Confidence | Runs | Agreement |
|---------------------|----------------|------------|------|-----------|
| crosslink_position  | start          | ✅ high    | 15   | 80.0%     |
| encode_eclip        | False          | ✅ high    | 15   | 80.0%     |
| skip_umi_dedupe     | False          | ✅ high    | 15   | 93.3%     |
| move_umi_to_header  | True           | ✅ high    | 15   | 73.3%     |
| umi_separator       | _              | ✅ high    | 15   | 73.3%     |
```

## Output Structure

```
clip_seq_output/
├── report.md              # Suggestion report with frequency tables
├── result.json            # Machine-readable suggestions + execution list
└── reproducibility/
    └── commands.sh        # Exact command to reproduce
```

## Gotchas

- **`crosslink_position` is protocol-specific, not universal.** iCLIP reads cross-link at the nucleotide immediately 5′ of the insert (use `start`); eCLIP reads have the cross-link at the 3′ end (use `end`). Mixing protocols in one history pool will dilute the signal — always filter by protocol (iCLIP vs eCLIP) before trusting the suggested value.

- **The `params` field is only populated in execution *detail* responses, not list responses.** Calling `/executions/owned` returns summary dicts with `params=None`. You must call `GET /executions/{id}` individually for each run to get the actual parameter values. The skill does this automatically, but it means analysing 50 runs requires 50 extra API calls.

- **New accounts have zero owned executions.** The skill supplements with public CLIP-Seq executions from flow.bio by default (`--no-public` to disable). Public run parameters reflect the broader community, which may differ from your lab's protocols — check `agreement_pct` and treat low-agreement suggestions with caution.

- **Boolean params are stored as strings `"true"`/`"false"` in the API**, not Python bools. The `extract_params()` function coerces them — but if you process the raw `result.json` externally, remember to parse them manually.

- **`paraclu_min_value` is sparse.** Most runs leave it unset (pipeline default). A `low` confidence suggestion here usually means one experimenter tried a custom value; it doesn't represent consensus.

## Safety

- **Local-first**: No genomic data is uploaded — only pipeline metadata (parameters, run status) is retrieved from the flow.bio API.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **No hallucinated parameters**: Every suggestion links back to an observed execution count. The skill will not invent parameter values not seen in the history.
- **Credentials**: Never logged to disk. The `commands.sh` reproducibility file redacts the password with `***`.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) fetches and aggregates.
The agent must NOT override the frequency-based suggestions with its own prior knowledge of CLIP-seq parameters.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- Query mentions "CLIP-seq", "iCLIP", "eCLIP", "PAR-CLIP" alongside "parameters", "settings", or "history"
- User asks "what should I set for my next CLIP run?"

**Chaining partners**:
- `flow-bio`: Launch the pipeline with the suggested parameters after this skill produces recommendations.
- `seq-wrangler`: Pre-QC the FASTQ before running the pipeline.
- `multiqc-reporter`: Aggregate QC results after the pipeline completes.

## Maintenance

- **Review cadence**: When the flow.bio CLIP-Seq pipeline releases a new major version (currently v1.7), check the parameter schema for new fields.
- **Staleness signals**: New pipeline parameters appear in `schema.inputs` at `/pipelines/versions/{id}`; update `PARAM_SCHEMA` in `clip_seq.py` accordingly.
- **Deprecation**: Archive to `skills/_deprecated/` if the flow.bio API changes its execution detail format.

## Citations

- goodwright/clipseq pipeline (flow.bio CLIP-Seq v1.7): https://github.com/goodwright/clipseq
- König et al. 2010, *iCLIP reveals the function of hnRNP particles in splicing at individual nucleotide resolution*, Nature Structural & Molecular Biology — original iCLIP protocol
- Van Nostrand et al. 2016, *Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP)*, Nature Methods — eCLIP protocol
- Smith et al. 2017, *UMI-tools: modeling sequencing errors in Unique Molecular Identifiers*, Genome Research — UMI deduplication
- Dobin et al. 2013, *STAR: ultrafast universal RNA-seq aligner*, Bioinformatics — STAR alignment
