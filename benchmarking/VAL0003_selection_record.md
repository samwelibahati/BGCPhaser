# VAL0003 prospective validation selection record

Status: terminal G5 A1.1 failure; excluded from G6-G9 and primary ranking-performance eligibility

Date frozen: 2026-08-27

## Benchmark identity

- Benchmark ID: `VAL0003`
- MIBiG accession: `BGC0000437`
- Organism: *Pseudomonas syringae* pv. *syringae* B728a
- NCBI Taxonomy ID: `205918`
- Biosynthetic class: `NRP`
- Reference accession: `CP000075.1`
- Assembly accession: `GCA_000012245.1`

## Frozen MIBiG-derived locus metadata

The exact row was taken from the locally generated `benchmarking/mibig_modular_pool.tsv` before download, assembly, G5 localization, candidate enumeration, or scoring.

- start coordinate: `3060225`
- end coordinate: `3182922`
- coordinate convention: MIBiG 1-based inclusive start / 1-based exclusive end
- locus span: `122697 bp`
- locus completeness: `Unknown`
- total module count: `9`

The arithmetic check `3182922 - 3060225 = 122697` is consistent with the frozen coordinate convention. These coordinates were not changed in response to the VAL0003 outcome.

## Prospective run selection

- selected run: `ERR005143`
- BioSample/sample accession: `SAMEA882553`
- secondary sample accession: `ERS000165`
- BioProject/study accession: `PRJEB2007`
- experiment accession: `ERX000536`
- scientific name: *Pseudomonas syringae* pv. *syringae* B728a
- library name: `PssB728a`
- library strategy: `WGS`
- library source: `GENOMIC`
- library selection: `RANDOM`
- library layout: `PAIRED`
- platform: `ILLUMINA`
- instrument: `Illumina Genome Analyzer`
- first public: `2010-02-26`
- candidate runs for BGC0000437 after frozen screening: `1`

Selection basis: `EXACT_STRAIN_NAME_AND_TAXID;STRAIN_SPECIFIC_GENOME_PROJECT;SINGLE_QUALIFYING_RUN;RANDOM_PAIRED_ILLUMINA_WGS`.

Selection was locked before assembly, candidate enumeration, intrinsic scoring, chemistry evaluation, truth ranking, or reference-similarity ranking.

Selection-lock hashes:

- prospective candidate table SHA-256: `35c35dad712f8833f50c2a476417fd2a24288c26b3f151718d85eae282d63c81`
- locked ENA one-row manifest SHA-256: `1ff93be34c3934bb935d65772884a2c91ff1e3a3eb76f9c9a3bb24220153a466`

## Verified raw FASTQs

Both ENA files completed on the first HTTPS attempt under downloader protocol `1.3-resumable-curl`. Observed byte counts and MD5 values exactly matched the predeclared ENA metadata before G4.

- `ERR005143_1.fastq.gz`: 130,401,291 bytes; ENA/observed MD5 `45bbba69aa6831d0754185596fb1b102`; SHA-256 `7cf3898805e5ed4ef91792bd4428f115c9822c52806b2980922d2adb730a53b2`
- `ERR005143_2.fastq.gz`: 138,851,038 bytes; ENA/observed MD5 `1217d44345bfc5104e869e0eee7a7003`; SHA-256 `9667e8b44fa94291f5df3bf2db88a8890dab03db5a389d231ccbe8716ff233a4`

Total compressed download size: `269252329 bytes`.

Download executable recorded by the provenance manifest: `/usr/bin/curl`, curl `8.7.1 (x86_64-apple-darwin25.0)`.

## G4 outcome

Frozen software:

- fastp `1.3.6`
- SPAdes `4.3.0`

G4 classification: `PASS_WITH_QC_WARNING`.

Integrity checks:

- raw SHA-256 manifest: `2/2 PASS`
- processed SHA-256 manifest: `4/4 PASS`
- canonical SPAdes SHA-256 manifest: `7/7 PASS`
- SHA-256 manifest problems: `0`
- `g4_integrity_status`: `PASS`

fastp summary:

- reads before filtering: `7102266`
- reads after filtering: `2196916`
- read retention fraction: `0.309326065794`
- bases before filtering: `255681576`
- bases after filtering: `76760170`
- base retention fraction: `0.300217838144`
- post-filter Q20 rate: `0.94392`
- post-filter Q30 rate: `0.866672`
- post-filter GC content: `0.560721`
- low-quality reads: `4986`
- reads with too many N: `0`
- too-short reads: `4900298`

Assembly summary:

- contigs: `5319`
- contig total bases: `3725011`
- contig maximum length: `14656`
- contig N50: `1129`
- scaffolds: `2437`
- scaffold total bases: `4061894`
- scaffold maximum length: `30845`
- scaffold N50: `3167`
- canonical GFA segments: `6833`
- canonical GFA links: `32`
- canonical GFA segment bases: `3740594`

SPAdes produced eleven warning lines. In addition to the macOS `setrlimit(2)` memory-limit warning, the log reported an improperly determined erroneous-connection coverage threshold and a k-mer coverage-model valley reset to 4. SPAdes nevertheless finished and produced every required canonical output. The approximately 30% read/base retention and assembly warnings are retained as QC observations. The frozen protocol contains no predeclared exclusion threshold for these observations, so no reassembly, parameter tuning, alternate k-mers, or rescue workflow was introduced.

## Pre-G5 implementation safeguard

Before any VAL0003 endpoint localization was inspected, the G5 implementation was hardened in two protocol-consistency respects:

1. if the unique left anchor ends and the unique right anchor begins on the same reference-forward oriented graph state, the left alignment must end at or before the right alignment begins on that state; otherwise G5 fails as coordinate-order incompatible;
2. a primary-validation G5 pass is recorded as `PENDING_LATER_GATES`, since final primary-validation eligibility still depends on G6-G9.

These safeguards did not alter the frozen 250-nt anchor length, >=95% query coverage threshold, >=98% identity threshold, minimap2 2.31 `asm5` preset, graph-walk construction, or physical-localization uniqueness rule. VAL0003 G5 remained uninspected when the implementation and regression tests were updated.

Focused regression test after the hardening:

- command: `python -m pytest tests/test_g5_contiguous_graph_walks.py -q`
- result: `4 passed in 0.23s`
- G5 implementation commit present before VAL0003 G5: `32213a3cc1310b5cb781f0fdfa05eab04c9d6de3`

## Exact reference frozen before G5

The reference FASTA was acquired from NCBI using accession `CP000075.1` and validated before any VAL0003 G5 endpoint localization.

- FASTA header: `>CP000075.1 Pseudomonas syringae pv. syringae B728a, complete genome`
- sequence length: `6093698 bp`
- FASTA SHA-256: `36a643d87d415d0ec19e159a89b2bd548c3e5e7bde52f34af7dca8e22092152b`
- independent `shasum -a 256` result: exact match
- reference validation: `PASS`

The locally frozen MIBiG interval `3060225-3182922` was used unchanged for VAL0003 G5.

## G5 A1.1 terminal outcome

Frozen G5 inputs and software:

- reference SHA-256: `36a643d87d415d0ec19e159a89b2bd548c3e5e7bde52f34af7dca8e22092152b`
- GFA SHA-256: `6901c705c8a52407d3b950c9b7d33c30f28c3618c1916187df8ca8a9447043cc`
- `contigs.paths` SHA-256: `9bfdfb8b8b26257d289b71f8131fcc20863652d9a3ad614a7d9707b9f1cda3ca`
- `scaffolds.paths` SHA-256: `cba83e42eb3c3a3090a9e736903ac3ff712cedb79028c235fa35b073c2aabefb`
- Python: `3.11.16`
- NetworkX: `3.5`
- minimap2 executable: `/Users/rmaghembe/miniforge3_arm64/pkgs/minimap2-2.31-h6bd33b9_0/bin/minimap2`
- minimap2 SHA-256: `007777300e7f1dc464300135d1b16e5c47bef40ae19a1017c1504d33b2509142`
- minimap2 version: `2.31-r1302`
- minimap2 preset: `asm5`
- anchor length: `250 nt`
- minimum query coverage: `0.95`
- minimum identity: `0.98`

Reference-forward anchor intervals used by the implementation:

- left: Python0 `3060224:3060474`
- right: Python0 `3182671:3182921`

A1.1 topology audit:

- assembler path records: `27324`
- missing-edge transition observations: `0`
- unique missing-edge transitions: `0`
- unique contiguous directed graph-walk chunks: `13662`
- total spelled contiguous-walk bases: `7481056`

Anchor-localization result:

- left physical localizations: `0`
- right physical localizations: `0`
- same-inner-state coordinate-order check: `NOT_EVALUATED`
- G5 status: `FAIL`
- recorded failure reason: `left_anchor_physical_localizations=0`
- primary-validation contribution: `NO_G5_FAILURE`

Both endpoint anchors therefore fail the frozen G5 uniqueness/localization requirement. The system stops at G5 and cannot proceed to G6, G7 candidate enumeration, G8 truth-representation testing, G9 scoreability, or ranking-performance analysis.

No threshold, anchor length, minimap2 preset, graph definition, reassembly strategy, or alternate reference was changed after observing this outcome. No rescue was attempted.

## Final outcome-masking state

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- candidate count: `NA`
- G6-G9: `NOT RUN`

VAL0003 remains in the complete prospective screening denominator and is reported as a terminal G5 candidate-localization failure. It does not contribute to the primary eligible ambiguous ranking-performance cohort.
