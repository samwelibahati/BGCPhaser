# VAL0002 prospective validation selection record

Status: frozen before G5 endpoint localization

Date frozen: 2026-08-26

## Benchmark identity

- Benchmark ID: `VAL0002`
- MIBiG accession: `BGC0000392`
- Organism: *Actinosynnema mirum* DSM 43827
- NCBI Taxonomy ID: `446462`
- Biosynthetic class: `NRP`
- Compound: mirubactin
- Reference accession: `CP001630.1`
- Assembly accession: `GCA_000023245.1`

## Frozen MIBiG-derived locus metadata

The exact row was taken from the locally generated `benchmarking/mibig_modular_pool.tsv` before G5 localization.

- start coordinate: `3113324`
- end coordinate: `3141040`
- coordinate convention: MIBiG 1-based inclusive start / 1-based exclusive end
- locus span: `27716 bp`
- locus completeness: `Unknown`
- NRPS module count: `3`
- PKS module count: `0`
- total modular count: `3`

The arithmetic check `3141040 - 3113324 = 27716` is consistent with the frozen coordinate convention. These coordinates, anchor length, localization thresholds, and mapper settings cannot be changed in response to the VAL0002 G5 result.

## Prospective run selection

- selected run: `SRR8526193`
- BioSample: `SAMN00001904`
- secondary sample: `SRS002071`
- BioProject: `PRJNA19705`
- experiment: `SRX5329028`
- library strategy: `WGS`
- library source: `GENOMIC`
- library selection: `RANDOM`
- library layout: `PAIRED`
- platform: `ILLUMINA`
- instrument: `Illumina HiSeq 2000`
- alternative same-sample run: `SRR3924417`
- alternative exclusion from selection: `library_selection=size fractionation`

Selection was locked before assembly, candidate enumeration, intrinsic scoring, chemistry evaluation, truth ranking, or reference-similarity ranking.

Selection-lock hashes:

- prospective candidate table SHA-256: `096a03c7df8791d3293987b17f2038ba0c12e60eb7ce741273acf02adda34d77`
- locked ENA one-row manifest SHA-256: `d14225d18396dd0281b915ef2654775f48ebab8398ebd40974ed7c690091e3ed`

## Verified raw FASTQs

Both mates matched the ENA byte counts and MD5 values before G4.

- `SRR8526193_1.fastq.gz`: 4,223,760,177 bytes; MD5 `3da274491d124a2322ddf24c3eaa2e4b`; SHA-256 `be2199247e994f0286c48a7cce545ba66eb58dfedabc626cdcd7001a836a66aa`
- `SRR8526193_2.fastq.gz`: 4,234,348,169 bytes; MD5 `a97ebe7fb29d12553b8ee2b46063609b`; SHA-256 `bf343b9b808d7de5522243f2b8493376785329e1ef8b2f1418ecdeef6ecc922a`

Download implementation: `1.3-resumable-curl`.

## G4 outcome

Frozen preprocessing and assembly software:

- fastp `1.3.6`
- SPAdes `4.3.0`

G4 integrity status: `PASS_WITH_QC_WARNING`.

Key recorded outputs:

- reads before filtering: `73119724`
- reads after filtering: `65726094`
- read retention fraction: `0.898883234297`
- bases before filtering: `10967958600`
- bases after filtering: `7325541156`
- base retention fraction: `0.667903793510`
- contigs: `22586`
- contig total bases: `16632354`
- contig N50: `2039`
- scaffolds: `22536`
- scaffold total bases: `16633565`
- scaffold N50: `2086`
- GFA segments: `23028`
- GFA links: `723`
- SHA-256 manifest problems: `0`

The eight captured warning lines are repeated macOS `setrlimit(2)` memory-limit warnings plus the SPAdes warning summary. SPAdes completed and produced every required canonical output. The fragmented/high-total-length assembly remains unchanged; no system-specific reassembly, subsampling, alternate k-mer series, or rescue workflow is permitted.

## Pre-G5 implementation freeze

The A1.1 contiguous directed graph-walk implementation was generalized to require an explicit analysis role before VAL0002 G5. The focused regression suite was run locally with Python 3.11 and pytest 8.4.2 before endpoint localization:

- `tests/test_g5_contiguous_graph_walks.py`: `3 passed`

VAL0002 must be invoked with `--analysis-role PRIMARY_VALIDATION`.

## Outcome masking state at this freeze

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- G5 endpoint localization inspected: `NO`

The next permitted operation is acquisition and checksum recording of the exact `CP001630.1` reference FASTA, followed by G5 using the frozen locus coordinates above.
