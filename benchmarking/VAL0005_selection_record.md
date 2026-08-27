# VAL0005 prospective validation selection record

Status: terminal G5 A1.1 failure; stopped before G6

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

The frozen interval was captured before read download, assembly, G5 localization, candidate enumeration, or scoring and was not altered in response to the G5 outcome.

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

Mate 1 required resumable HTTPS retries and verified on attempt 16; mate 2 verified on attempt 1. The transport history is retained in `download_diagnostics.tsv` and did not alter finalized FASTQ integrity.

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

SPAdes emitted ten warning-containing log lines, all attributable to the repeated macOS `setrlimit(2)` failure to impose the requested 250-Gb process memory limit plus the pipeline summary that warnings were present. No coverage-model, erroneous-connection-threshold, graph-simplification, or assembly-content warning was reported. SPAdes completed and all canonical outputs passed SHA-256 verification. No reassembly or parameter change was introduced.

## Exact reference frozen before G5

The exact reference FASTA was acquired by accession `FN667742.1` and validated before any VAL0005 endpoint localization.

- FASTA header: `>FN667742.1 Xenorhabdus nematophila ATCC 19061 chromosome, complete genome`
- sequence length: `4432590 bp`
- FASTA SHA-256: `624174a6d371828bf52fd11a264a90835054d7c67914f2f844d16da840521a45`
- independent `shasum -a 256` result: exact match
- reference validation: `PASS`

## G5 A1.1 outcome

G5 used the exact frozen reference and locus interval with minimap2 `2.31-r1302`, preset `asm5`, 250-nt endpoint anchors, minimum query coverage `0.95`, and minimum identity `0.98`.

Canonical G5 input hashes:

- reference FASTA SHA-256: `624174a6d371828bf52fd11a264a90835054d7c67914f2f844d16da840521a45`
- GFA SHA-256: `37bb9daebd737a6d2f270cdd4c8deca1f7b92a9c0d8085bfecc348962aaf6b3e`
- `contigs.paths` SHA-256: `0b84ceeda4f732e4148d668106702e83ee45ccaa213b97cb0541ac9783af0acc`
- `scaffolds.paths` SHA-256: `6ffcbc004a38b8a834bf6b3077422035606d7181457c9d181398d7acaacf2d7c`
- minimap2 executable SHA-256: `007777300e7f1dc464300135d1b16e5c47bef40ae19a1017c1504d33b2509142`

Graph-walk localization summary:

- assembler path records: `2612`
- missing-edge transition observations: `0`
- unique missing-edge transitions: `0`
- unique contiguous directed graph-walk chunks: `1306`
- total spelled contiguous-walk bases: `8953488`
- left physical localization count: `0`
- right physical localization count: `0`
- same-inner-state coordinate-order test: `NOT_EVALUATED`

G5 classification: `FAIL` / `G5_A1_1_FAIL`.

Both predefined 250-nt endpoint anchors lacked a qualifying physical localization under the frozen A1.1 rule. The implementation reports the first terminal reason as `left_anchor_physical_localizations=0`; the summary independently records `right_physical_localization_count=0` as well. No G6 boundary states were emitted. Primary-validation contribution is `NO_G5_FAILURE`.

VAL0005 therefore stops at G5 and does not proceed to G6, candidate enumeration, intrinsic scoring, chemistry evaluation, truth ranking, or reference-similarity analysis. No threshold, anchor length, minimap2 preset, reference, assembly, or rescue rule was altered after observing this outcome.

## Outcome masking state at terminal stop

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- G5 endpoint localization inspected: `YES — terminal gate outcome only`

No further primary-validation analysis is permitted for VAL0005 beyond preservation of the frozen audit artifacts.
