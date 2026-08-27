# VAL0005 prospective validation selection record

Status: G4 passed with environmental QC warning; exact reference verified; G5 not yet run

Date frozen: 2026-08-27

## Benchmark identity

- Benchmark ID: `VAL0005`
- MIBiG accession: `BGC0000416`
- Organism: *Xenorhabdus nematophila* ATCC 19061
- NCBI Taxonomy ID: `406817`
- Biosynthetic class: `NRP`
- Reference accession: `FN667742.1`
- Assembly accession: `GCA_000252955.1`

## Frozen MIBiG-derived locus metadata

- start coordinate: `2154724`
- end coordinate: `2170060`
- coordinate convention: MIBiG 1-based inclusive start / 1-based exclusive end
- locus span: `15336 bp`
- locus completeness: `complete`
- total module count: `4`

The frozen interval was captured before read download, assembly, G5 localization, candidate enumeration, or scoring and cannot be altered in response to later outcomes.

## Prospective run selection

- selected run: `SRR21665144`
- BioSample/sample accession: `SAMN16789554`
- secondary sample accession: `SRS15195604`
- BioProject/study accession: `PRJNA678184`
- experiment accession: `SRX17663552`
- scientific name: *Xenorhabdus nematophila* ATCC 19061
- sample title: `MIGS Cultured Bacterial/Archaeal sample from Xenorhabdus nematophila ATCC 19061`
- experiment title: `Illumina MiSeq sequencing: Illumina TruSeq PCR-free library of Xenorhabdus nematophila ATCC 19061: 2x 300nt paired-end run`
- library strategy: `WGS`
- library source: `GENOMIC`
- library selection: `RANDOM`
- library layout: `PAIRED`
- platform: `ILLUMINA`
- instrument: `Illumina MiSeq`
- candidate runs for BGC0000416 after frozen screening: `1`

Selection basis: `EXACT_TYPE_STRAIN_NAME_AND_TAXID;FINISHED_REFERENCE_SAME_ATCC_STRAIN;SINGLE_QUALIFYING_RUN;RANDOM_PAIRED_ILLUMINA_WGS`.

This selection is a defensibly matched same named ATCC strain. It does not assert that the 2026 sequencing BioSample is the same historical BioSample used to construct the finished reference.

Selection-lock hashes:

- prospective candidate table SHA-256: `35c35dad712f8833f50c2a476417fd2a24288c26b3f151718d85eae282d63c81`
- locked ENA one-row manifest SHA-256: `dd1357a46dec87b4933902526149e1dcff467316548d749987aab3e7dae3e183`

## Verified raw FASTQs

Both files were finalized only after exact ENA byte-count and MD5 verification under downloader protocol `1.3-resumable-curl` using `/usr/bin/curl` 8.7.1.

- `SRR21665144_1.fastq.gz`: 252,462,925 bytes; MD5 `53494d1b8e767d87f1a07288f295fadd`; SHA-256 `bfee42f320bbaa85f245c0aeeb36b810b609e91ee074e53a115bd144f1df0a92`
- `SRR21665144_2.fastq.gz`: 315,890,191 bytes; MD5 `f520582f699c11e7d2eba972445e2174`; SHA-256 `3e50269b3b08d4c7295d29d284b6241639ecd903aa566128f3db6a6bf4dc0514`
- total verified compressed bytes: `568353116`

Mate 1 required resumable HTTPS retries and verified on attempt 16; mate 2 verified on attempt 1. The transport history is retained in `download_diagnostics.tsv` and does not alter finalized FASTQ integrity.

## G4 outcome

Frozen software:

- fastp `1.3.6`
- SPAdes `4.3.0`

G4 classification: `PASS_WITH_QC_WARNING`.

All integrity manifests passed:

- raw SHA-256 manifest: `2/2 PASS`
- processed SHA-256 manifest: `4/4 PASS`
- canonical SPAdes SHA-256 manifest: `7/7 PASS`
- SHA-256 manifest problems: `0`
- `g4_integrity_status`: `PASS`

fastp summary:

- reads before filtering: `3060286`
- reads after filtering: `2954414`
- read retention fraction: `0.965404540621`
- bases before filtering: `899699637`
- bases after filtering: `791839584`
- base retention fraction: `0.880115486809`
- post-filter Q20 rate: `0.99267`
- post-filter Q30 rate: `0.959816`
- post-filter GC content: `0.43653`
- low-quality reads: `288`
- reads with too many N: `0`
- too-short reads: `105250`

Assembly summary:

- contigs: `653`
- contig total bases: `4475982`
- contig maximum length: `237637`
- contig N50: `60384`
- scaffolds: `652`
- scaffold total bases: `4476082`
- scaffold maximum length: `237637`
- scaffold N50: `60384`
- canonical GFA segments: `1735`
- canonical GFA links: `2332`
- canonical GFA segment bases: `4601711`
- GFA bytes: `4754738`
- `contigs.paths` bytes: `75741`
- `scaffolds.paths` bytes: `75670`

SPAdes emitted ten warning-containing log lines, all attributable to the repeated macOS `setrlimit(2)` failure to impose the requested 250-Gb process memory limit plus the pipeline summary that warnings were present. No coverage-model, erroneous-connection-threshold, graph-simplification, or assembly-content warning was reported. SPAdes completed and all canonical outputs passed SHA-256 verification. The warning is therefore retained transparently as an environmental QC observation; no reassembly, parameter change, or exclusion is introduced.

## Exact reference frozen before G5

The exact reference FASTA was acquired by accession `FN667742.1` and validated before any VAL0005 endpoint localization.

- FASTA header: `>FN667742.1 Xenorhabdus nematophila ATCC 19061 chromosome, complete genome`
- sequence length: `4432590 bp`
- FASTA SHA-256: `624174a6d371828bf52fd11a264a90835054d7c67914f2f844d16da840521a45`
- independent `shasum -a 256` result: exact match
- reference validation: `PASS`

The frozen MIBiG interval `2154724-2170060` lies within this exact sequence and remains the only permitted endpoint interval for VAL0005 G5.

## Outcome masking state at this freeze

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- G5 endpoint localization inspected: `NO`

The next permitted operation is G5 A1.1 anchor localization using this exact reference, the canonical G4 graph/path files, frozen interval `2154724-2170060`, minimap2 2.31 `asm5`, 250-nt endpoint anchors, >=95% query coverage and >=98% identity. No rescue rule, threshold change, reassembly, alternate reference, or outcome-driven substitution is permitted.
