# BGCPhaser external-validation assembly protocol v1

Status: frozen before the first new validation assembly

## Purpose

Generate a standardized short-read assembly graph for every prospectively qualified paired-end Illumina isolate dataset entering BGCPhaser expanded validation.

The same preprocessing and assembly settings are applied to every new validation system. Graph ambiguity, candidate count and BGCPhaser score are unavailable when these settings are applied.

## Software

- fastp v1.3.6
- SPAdes v4.3.0

Exact executable versions are written to each benchmark provenance record before processing.

## Input contract

Each system must provide two raw FASTQ files representing the complete paired-end Illumina run selected at gate G3. Raw-file MD5 or SHA-256 digests are recorded before preprocessing.

Runs containing only one mate, merged-only reads, amplicon libraries, targeted-enrichment libraries, metagenomes or synthetic/simulated reads do not enter this workflow.

## Read preprocessing

fastp is run once per selected run. The locked intent is adapter removal, Q20 right-end quality trimming and retention of surviving paired and singleton reads.

Command template:

```bash
fastp \
  --in1 RAW_R1.fastq.gz \
  --in2 RAW_R2.fastq.gz \
  --out1 clean_R1.fastq.gz \
  --out2 clean_R2.fastq.gz \
  --unpaired1 unpaired_R1.fastq.gz \
  --unpaired2 unpaired_R2.fastq.gz \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --length_required 30 \
  --thread 4 \
  --json fastp.json \
  --html fastp.html
```

No deduplication, read merging, base correction, fixed-length truncation or coverage subsampling is permitted in the primary validation workflow.

Both singleton files are retained for downstream direct read-junction evidence. They do not contribute to paired-fragment linkage evidence.

## Canonical SPAdes assembly

The canonical graph is generated from the surviving paired reads with SPAdes isolate mode:

```bash
spades.py \
  --isolate \
  -1 clean_R1.fastq.gz \
  -2 clean_R2.fastq.gz \
  -o spades
```

SPAdes automatic k-mer selection is retained. No reference genome, MIBiG sequence, long reads, trusted contigs, untrusted contigs, coverage cutoff override or manual k-mer series enters the primary assembly.

The primary assembly is invalid if SPAdes terminates unsuccessfully or does not produce the graph/path files required by BGCPhaser.

## Canonical outputs

The following SPAdes outputs are retained when present:

- `assembly_graph_with_scaffolds.gfa`
- `assembly_graph_after_simplification.gfa`
- `contigs.paths`
- `scaffolds.paths`
- `contigs.fasta`
- `scaffolds.fasta`
- `params.txt`
- `spades.log`

For BGCPhaser candidate construction the canonical graph is `assembly_graph_with_scaffolds.gfa`, accompanied by the frozen union of `contigs.paths` and `scaffolds.paths`.

## Provenance

For each benchmark record:

1. raw FASTQ digests;
2. processed FASTQ digests;
3. fastp version and command;
4. SPAdes version and command;
5. SHA-256 digests of the canonical graph, contigs/scaffolds and path files;
6. read counts before and after filtering;
7. total bases before and after filtering;
8. SPAdes exit status and warnings.

These are recorded before BGC anchor localization or candidate enumeration results are inspected.

## Failure rule

A system that cannot complete this fixed workflow fails gate G4 and remains in the benchmark registry with the failure reason. No system-specific trimming, alternate assembler, alternate k-mer series or rescue assembly is introduced for primary-validation eligibility.
