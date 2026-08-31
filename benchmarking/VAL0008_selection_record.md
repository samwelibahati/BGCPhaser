# VAL0008 prospective validation selection record

Status: G3 passed; raw paired-end reads downloaded and integrity-verified; G4 not yet run

Selection frozen: 2026-08-28
Raw-read verification completed: 2026-08-30

## Benchmark identity

- Benchmark ID: `VAL0008`
- MIBiG accession: `BGC0001479`
- Organism: *Nostoc punctiforme* PCC 73102
- NCBI Taxonomy ID: `63737`
- Biosynthetic class: `NRP`
- Reference accession: `CP001037.1`
- Assembly accession: `GCA_000020025.1`

## Frozen MIBiG-derived locus metadata

- start coordinate: `3052001`
- end coordinate: `3081156`
- locus span: `29155 bp`
- locus completeness: `complete`
- total module count: `6`

The locus and candidate ordering were frozen before read download, assembly, G5 localization, candidate enumeration, or BGCPhaser scoring.

## Prospective run selection

- selected run: `SRR30545090`
- BioSample/sample accession: `SAMN43372961`
- secondary sample accession: `SRS22462789`
- BioProject/study accession: `PRJNA1152567`
- experiment accession: `SRX25968490`
- scientific name: *Nostoc punctiforme* PCC 73102
- sample title: `N. punctiforme cells, Dark24, bio rep 1`
- experiment title: `NextSeq 500 sequencing: WGS of Nostoc punctiforme Dark24`
- library strategy: `WGS`
- library source: `GENOMIC`
- library selection: `PCR`
- explicitly expected library selection: `PCR`
- library layout: `PAIRED`
- platform: `ILLUMINA`
- instrument: `NextSeq 500`
- read count: `97604601`
- base count: `30157965721`
- provenance status: `NO_DERIVATIVE_KEYWORD_DETECTED`
- provenance flag: empty

Selection basis: `EXACT_STRAIN_NAME_AND_TAXID;COMPLETE_REFERENCE_SAME_PCC_STRAIN;PLASMID_FOCUSED_RUNS_EXCLUDED;WHOLE_GENOME_PAIRED_ILLUMINA_RUN_SELECTED;PCR_LIBRARY_SELECTION_EXPLICITLY_PREDECLARED`.

The six alternative PRJNA1034103 rows (`SRR26619834` through `SRR26619839`) were excluded prospectively because their sample/experiment metadata describe plasmid-focused sequence from desiccated *N. punctiforme* material. `SRR30545090` was selected before any downstream outcome was inspected because its experiment is explicitly whole-genome sequencing of PCC 73102. `library_selection=PCR` was prospectively accepted as an exact metadata value for this WGSeq run; it is not being reinterpreted as RANDOM.

## Lock implementation provenance

The first lock attempt failed before creating any lock output because `benchmarking/scripts/12_lock_validation_run.py` had historically hard-coded `library_selection=RANDOM`. This was an implementation constraint stricter than the frozen G3 protocol, which requires paired-end Illumina WGS and excludes amplicon/targeted-enrichment datasets but does not require the literal ENA value RANDOM.

Before rerunning the lock, the helper was changed narrowly so that `--expected-library-selection` is explicit and defaults to `RANDOM`, preserving all prior lock behavior. The change was committed before the successful VAL0008 lock:

- helper-change commit: `861c0e724ba1a79bc0775c28ee7bfb171f86b3fb`
- helper blob/content SHA after change: `05212e43a93eebae78c28a8aa960fdbde70a05de`

No biological eligibility criterion, score, threshold, reference, or downstream gate was modified.

## Frozen lock hashes

- prospective candidate table SHA-256: `35c35dad712f8833f50c2a476417fd2a24288c26b3f151718d85eae282d63c81`
- locked one-row ENA manifest SHA-256: `dd92e4a2b34a47086e24f4c6c2cbcaaeccae27cc733b11eae808f13990c18422`
- lock status: `PASS`

## Verified raw FASTQ provenance

The exact locked ENA FASTQs were acquired with downloader protocol `1.3-resumable-curl` using `/usr/bin/curl` 8.7.1. Transfers were interrupted repeatedly and resumed without deleting the retained partial files. Finalization occurred only after exact expected byte count and ENA MD5 verification.

### Mate 1

- file: `SRR30545090_1.fastq.gz`
- expected bytes: `9279646920`
- observed bytes: `9279646920`
- ENA MD5: `49a2910a3bb97e10ea0bee14193559e4`
- observed MD5: `49a2910a3bb97e10ea0bee14193559e4`
- SHA-256: `8ee45de3830e9b1161783d53a8ea0c83c6d774fdfd60c57b9458b67c8e32dc99`
- final provenance transport field: `existing_verified`

`existing_verified` records that the final successful downloader invocation found the already-completed mate-1 file and reverified it before processing mate 2; it does not imply an unverified pre-existing input.

### Mate 2

- file: `SRR30545090_2.fastq.gz`
- expected bytes: `9610960013`
- observed bytes: `9610960013`
- ENA MD5: `c98d509a8e667b44e55c40270671ee95`
- observed MD5: `c98d509a8e667b44e55c40270671ee95`
- SHA-256: `39ae61f47cfbb22aae91e091d68ee5a4b7e93b7c060781a56453d4958f451f02`
- final transport: `curl_ftp_resumable`

### Download-audit preservation

The terminal verification showed the final mate-2 FTP attempt reaching the exact expected byte count with curl exit code 0 and result `VERIFIED`. The raw directory contains:

- `download_diagnostics.tsv`
- `download_diagnostics_pre_resume.tsv`
- `download_diagnostics_pre_resume_2.tsv`
- `download_diagnostics_pre_resume_3.tsv`

No `*.part` files remained after verification. Thus the complete compressed raw input set is exactly `18890606933` bytes and both raw FASTQs are integrity-verified before G4.

## Outcome masking state after raw-read verification

- BGCPhaser score inspected: `NO`
- truth rank inspected: `NO`
- chemistry outcome inspected: `NO`
- reference-similarity outcome inspected: `NO`
- G5 endpoint localization inspected: `NO`
- G4 assembly outcome inspected: `NO`

The next permitted operation is the frozen G4 fastp/SPAdes workflow. No alternate run, read subset, reference substitution, preprocessing change, or downstream-outcome-driven modification is permitted.
