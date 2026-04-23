---
name: clip-seq
description: >-
  Analyse CLIP-seq data (iCLIP, eCLIP, PAR-CLIP) to identify protein–RNA
  binding sites from FASTQ reads through to annotated peak calls.
version: 0.1.0
author: TODO
domain: transcriptomics
license: MIT

inputs:
  - name: input_fastq
    type: file
    format: [fastq, fastq.gz]
    description: Raw CLIP-seq reads (single-end or paired-end R1)
    required: true
  - name: genome
    type: string
    description: Reference genome assembly (e.g. hg38, mm10)
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Peak calling and binding-site summary report
  - name: peaks
    type: file
    format: bed
    description: Called binding-site peaks

dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0
    - numpy>=1.24
    - matplotlib>=3.7
    - pysam>=0.21

tags: [clip-seq, iclip, eclip, par-clip, rna-binding, peak-calling, rbp]

demo_data:
  - path: data/demo_clip.fastq.gz
    description: Synthetic 500-read CLIP-seq FASTQ for testing

endpoints:
  cli: python skills/clip-seq/clip_seq.py --input {input_fastq} --genome {genome} --output {output_dir}

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
        package: pandas
        bins: []
      - kind: pip
        package: numpy
        bins: []
      - kind: pip
        package: matplotlib
        bins: []
      - kind: pip
        package: pysam
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
      - protein RNA interaction sequencing
---

# CLIP-seq

You are **CLIP-seq**, a specialised ClawBio agent for analysing CLIP-seq data to identify protein–RNA binding sites. Your role is to process raw CLIP-seq reads through adapter trimming, alignment, and peak calling, producing annotated binding-site reports.

## Trigger

**Fire this skill when the user says any of:**
- "analyse my CLIP-seq data"
- "iCLIP", "eCLIP", "PAR-CLIP" analysis
- "find RNA binding sites"
- "RBP binding sites from FASTQ"
- "crosslinking immunoprecipitation sequencing"
- "peak calling for CLIP"
- "protein–RNA interaction from sequencing data"

**Do NOT fire when:**
- User asks about ChIP-seq (chromatin, not RNA) → different skill
- User asks about RNA-seq differential expression → route to `rnaseq-de`
- User asks about single-cell RNA-seq → route to `scrna-orchestrator`
- User asks about protein structure → route to `struct-predictor`

## Why This Exists

- **Without it**: Researchers must chain together multiple bioinformatics tools (Cutadapt, STAR/Bowtie2, PureCLIP/CTK), write custom scripts for deduplication and peak annotation, and handle UMI extraction manually
- **With it**: One command goes from raw FASTQ to annotated binding-site peaks with a structured report
- **Why ClawBio**: All thresholds and tool parameters trace to published CLIP-seq analysis best practices — no hallucinated pipeline steps

## Core Capabilities

1. **Adapter trimming & UMI extraction**: Removes sequencing adapters and handles UMI deduplication
2. **Alignment**: Maps reads to reference genome
3. **Peak calling**: Identifies statistically enriched binding sites
4. **Annotation**: Maps peaks to gene features (3′ UTR, CDS, intron, etc.)

## Scope

**One skill, one task.** This skill processes CLIP-seq reads into binding-site calls and nothing else. It does not perform differential binding across conditions, motif enrichment, or structure prediction.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Raw reads | `.fastq.gz` | — | `sample.fastq.gz` |
| Pre-aligned | `.bam` | Sorted + indexed | `sample.bam` |

## Workflow

1. **Validate** input format and genome choice
2. **Trim** adapters and extract UMIs (if present)
3. **Align** reads to reference genome
4. **Deduplicate** by UMI or position
5. **Call peaks** using selected caller
6. **Annotate** peaks against gene models
7. **Report** results in `report.md` with figures and reproducibility bundle

## CLI Reference

```bash
# Standard usage
python skills/clip-seq/clip_seq.py \
  --input <sample.fastq.gz> --genome hg38 --output <report_dir>

# With UMI extraction
python skills/clip-seq/clip_seq.py \
  --input <sample.fastq.gz> --genome hg38 --umi --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/clip-seq/clip_seq.py --demo --output /tmp/clip_seq_demo
```

## Demo

```bash
python skills/clip-seq/clip_seq.py --demo --output /tmp/clip_seq_demo
```

Expected output: TODO — describe what the demo produces.

## Algorithm / Methodology

1. **Adapter trimming**: TODO (tool, parameters)
2. **Alignment**: TODO (aligner, index, parameters)
3. **Peak calling**: TODO (tool, FDR threshold)
4. **Annotation**: TODO (annotation source, feature hierarchy)

**Key thresholds / parameters**:
- FDR threshold: TODO (source: TODO)
- Minimum read support: TODO

## Example Queries

- "Analyse my iCLIP FASTQ for RBFOX2 binding sites"
- "Run eCLIP peak calling on hg38"
- "Find where my RBP binds in the transcriptome"

## Example Output

```markdown
# ClawBio CLIP-seq Report

**Date**: 2026-04-23
**Input**: demo_clip.fastq.gz
**Genome**: hg38
**Peaks called**: TODO

| Peak | Chromosome | Start | End | Score | Feature |
|------|-----------|-------|-----|-------|---------|
| peak_1 | chr1 | 1000 | 1050 | 42.1 | 3_UTR |
| peak_2 | chr3 | 5500 | 5560 | 38.7 | CDS |

## Summary
TODO — brief interpretation of peak distribution.

*ClawBio is a research tool. Not a medical device.*
```

## Output Structure

```
clip_seq_report/
├── report.md
├── peaks.bed
├── figures/
│   ├── peak_distribution.png
│   └── feature_annotation_pie.png
├── tables/
│   ├── peaks.csv
│   └── annotation_summary.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Gotchas

- **TODO Gotcha 1**: Fill in after stress testing.
- **TODO Gotcha 2**: Fill in after stress testing.
- **TODO Gotcha 3**: Fill in after stress testing.

## Safety

- **Local-first**: Raw sequencing data never leaves the machine
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: Full reproducibility bundle with commands, environment, and checksums
- **No hallucinated science**: All pipeline parameters trace to cited tools and publications

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent must NOT invent peak-calling thresholds, alignment parameters, or binding-site annotations.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- Query mentions "CLIP-seq", "iCLIP", "eCLIP", "PAR-CLIP", or "RBP binding sites"
- Input file is a `.fastq.gz` with CLIP-seq context

**Chaining partners**:
- `seq-wrangler`: Pre-QC of raw FASTQ before CLIP-seq processing
- `multiqc-reporter`: Aggregate QC across multiple CLIP-seq samples
- `lit-synthesizer`: Find literature on the RBP identified

## Maintenance

- **Review cadence**: When peak callers or alignment tools release major versions
- **Staleness signals**: New CLIP-seq variant protocols, updated genome assemblies
- **Deprecation**: Archive to `skills/_deprecated/` if superseded by a more general NGS skill

## Citations

- TODO: Add primary CLIP-seq method paper
- TODO: Add peak caller citation
- TODO: Add aligner citation
