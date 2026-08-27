# VAL0003 prospective validation selection record

Status: G3 passed; raw paired-end reads verified; G4 not yet run

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

The arithmetic check `3182922 - 3060225 = 122697` is consistent with the frozen coordinate convention. These coordinates cannot be changed in response to later VAL0003 outcomes.

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

## Outcome masking state at this freeze

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- G5 endpoint localization inspected: `NO`

The next permitted operation is frozen G4 preprocessing and SPAdes assembly using fastp 1.3.6 and SPAdes 4.3.0. No system-specific tuning, subsampling, reference guidance, trusted contigs, or alternative rescue assembly is permitted.
