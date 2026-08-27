# VAL0005 prospective validation selection record

Status: G3 passed; raw paired-end reads downloaded and verified; G4 not yet run

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

## Predeclared FASTQ files

- `SRR21665144_1.fastq.gz`: expected 252,462,925 bytes; ENA MD5 `53494d1b8e767d87f1a07288f295fadd`
- `SRR21665144_2.fastq.gz`: expected 315,890,191 bytes; ENA MD5 `f520582f699c11e7d2eba972445e2174`
- total expected compressed download: `568353116 bytes`

## Verified raw FASTQs

Both files were finalized only after exact ENA byte-count and MD5 verification under downloader protocol `1.3-resumable-curl` using `/usr/bin/curl` 8.7.1.

- `SRR21665144_1.fastq.gz`: observed 252,462,925 bytes; observed MD5 `53494d1b8e767d87f1a07288f295fadd`; SHA-256 `bfee42f320bbaa85f245c0aeeb36b810b609e91ee074e53a115bd144f1df0a92`
- `SRR21665144_2.fastq.gz`: observed 315,890,191 bytes; observed MD5 `f520582f699c11e7d2eba972445e2174`; SHA-256 `3e50269b3b08d4c7295d29d284b6241639ecd903aa566128f3db6a6bf4dc0514`
- total verified compressed bytes: `568353116`

Mate 1 encountered transient HTTPS transport interruptions. Attempts 1-6 ended with curl exit code 52 before bytes were retained; attempt 7 reached 240,074,752 bytes before curl exit code 18; attempts 8-13 retained that partial file with exit code 52; attempts 14 and 15 advanced the partial file to 240,138,435 and 240,210,118 bytes respectively with exit code 18; attempt 16 resumed from 240,210,118 bytes and completed with the exact expected byte count and MD5. Mate 2 completed and verified on HTTPS attempt 1. These transport interruptions did not alter the finalized FASTQ contents and are preserved in `download_diagnostics.tsv` as provenance rather than treated as biological or sequencing QC failures.

## Outcome masking state at this freeze

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- G5 endpoint localization inspected: `NO`

The next permitted operation is frozen G4 preprocessing and SPAdes assembly using fastp 1.3.6 and SPAdes 4.3.0. No alternate run, read subset, assembly parameter tuning, reference guidance, trusted contigs, or outcome-driven substitution is permitted.
