# Remaining prospective candidate screen — 2026-08-27

Status: frozen before selection/download of the next post-VAL0005 validation system

## Purpose

This record freezes the locus-level review of the remaining high-ranked rows in `benchmarking/val0003_provenance_candidates.tsv` after terminal G5 failure of VAL0005. The review was performed without inspecting any BGCPhaser candidate score, truth rank, chemistry outcome, reference-similarity outcome, or G5 localization outcome for any unprocessed system.

Prospective candidate-table SHA-256 already frozen by prior locks: `35c35dad712f8833f50c2a476417fd2a24288c26b3f151718d85eae282d63c81`.

The frozen validation source universe is bacterial modular NRPS, type I PKS, or hybrid NRPS-PKS loci with complete/high-confidence finished sequence-level truth and real paired-end Illumina WGS from the same isolate or a predeclared defensibly matched biological source. Screening stops at the first failed gate.

## Frozen decisions in ranked order

### BGC0002060 — duplicate biological system; no new independent benchmark ID

- Organism: *Pseudomonas syringae* pv. *syringae* B728a
- Reference/run: `CP000075.1` / `ERR005143`
- Same biological source, reference and raw run as VAL0003.
- Decision: retain only as a within-genome secondary locus; it cannot contribute an additional independent biological system.

### BGC0000336 — same NRRL 11379 biological system as VAL0004; G2-incompatible locus truth

- Organism: *Streptomyces filamentosus* NRRL 11379
- Reference: `AY787762.1`
- The accession is a partial daptomycin biosynthetic-cluster sequence, not a complete/high-confidence finished genome reference.
- The four candidate runs are the same NRRL 11379 sequencing project already encountered for VAL0004.
- Decision: no new independent benchmark; the locus is independently incompatible with G2.

### BGC0001054 / BGC0001132 / BGC0001692 — duplicate biological system; no new independent benchmark IDs

- Organism: *Xenorhabdus nematophila* ATCC 19061
- Run: `SRR21665144`
- Reference chromosome represented by `FN667742.1` or its RefSeq counterpart `NC_014228.1`.
- Same biological source and raw run as VAL0005.
- Decision: retain only as within-genome secondary loci; they cannot inflate the independent-system count.

### BGC0001290 / BGC0001399 / BGC0001515 / BGC0001699 — outside frozen source universe

- Organism: *Aspergillus nidulans* FGSC A4.
- Decision: excluded before prospective bacterial benchmark assignment. The frozen source universe is explicitly bacterial; no fungal locus is promoted without a separate prospective protocol expansion.
- The ISS/ground growth metadata are therefore not used to select among these fungal rows.

### BGC0001443 — distinct bacterial system; G2 FAIL

- Organism: *Actinomadura atramentaria* DSM 43919
- MIBiG reference: `NZ_KB907224.1`
- Candidate run: `SRR896000`; exact strain, WGS/GENOMIC/RANDOM/PAIRED/ILLUMINA.
- The reference belongs to assembly ASM38188v1 / `GCA_000381885` and the genome product is explicitly a scaffold-level Microbial Minimal Draft Isolate Genome.
- Decision: G2 fails before run locking or read download. MIBiG locus completeness does not convert a draft/scaffold genome into finished sequence-level truth.
- Benchmark identifier reserved: `VAL0006`.

### BGC0001467 — distinct bacterial system; G2 FAIL

- Organism: *Fischerella* sp. PCC 9431
- MIBiG reference: `NZ_KE650771.1`
- Two candidate Illumina WGS runs were present: `SRR3937996` and `SRR610343`.
- The reference is part of WGS master `ALVX00000000.1`; the public assembly contains 43 contigs organized into three scaffolds and is described as an improved high-quality draft.
- Decision: G2 fails before adjudicating the two-run choice or downloading reads.
- Benchmark identifier reserved: `VAL0007`.

### BGC0001479 — next distinct system passing the present pre-download screen

- Organism: *Nostoc punctiforme* PCC 73102, taxid `63737`
- Frozen MIBiG reference: `CP001037.1`
- Frozen local locus interval: `3052001-3081156`
- Locus class: NRPS
- Locus completeness: complete
- Modules: 6
- `CP001037.1` is the complete PCC 73102 chromosome; assembly `GCA_000020025.1` is a complete genome.
- Seven candidate run rows were present.
- `SRR26619834` through `SRR26619839` belong to PRJNA1034103 (`Nostoc space cargo`) and their metadata explicitly describe plasmid sequence from extracted/desiccated *N. punctiforme*. The associated publication states that PRJNA1034103 contains the plasmid sequences used in that study. These six rows fail G3 as targeted/plasmid-derived data and are not eligible alternatives.
- `SRR30545090`, sample `SAMN43372961`, study `PRJNA1152567`, is labeled WGS/GENOMIC/PAIRED/ILLUMINA NextSeq 500 for *N. punctiforme* PCC 73102. The associated study explicitly identifies `SRR30545090` as whole-genome-sequencing (`WGSeq`) raw reads. Its `library_selection=PCR` field is not treated as targeted enrichment: the experiment is independently documented as whole-genome sequencing and the frozen G3 exclusions concern amplicon/targeted-enrichment data, not generic PCR amplification of a WGS library.
- The biological source is PCC 73102 / ATCC 29133; NCBI recognizes these as co-identical strain designations.
- Decision: `SRR30545090` is the predeclared next run to attempt G3 lock. The six plasmid-focused runs are excluded before lock. No other run substitution is permitted in response to later assembly or G5 outcomes.
- Next benchmark identifier if the lock passes: `VAL0008`.

### BGC0001692 — duplicate of VAL0005

Already covered above with the ATCC 19061 within-genome loci.

### BGC0001796 — eligible-looking reserve, not selected ahead of BGC0001479

- Organism: *Brevibacillus laterosporus* DSM 25
- Reference: `CP017705.1`, publicly designated complete genome
- Candidate run: `SRR892139`; exact strain, WGS/GENOMIC/RANDOM/PAIRED/ILLUMINA HiSeq 2000; one candidate row.
- Decision: retain as the next high-priority independent reserve after the BGC0001479 system. It is not advanced now so that ranked screening order is preserved.

### BGC0001801 — G2-incompatible frozen truth

- Organism: *Streptomyces venezuelae* ATCC 10712
- Frozen reference: `FR845719.1`
- Although the accession is labeled a complete genome, later completion work documented 9,156 undetermined sequence characters in FR845719 and reported that the replacement complete sequence `CP029197` filled 147 gaps and corrected 46 overlapping regions.
- Decision: the frozen MIBiG truth accession does not meet the high-confidence finished-reference standard. G2 fails. `CP029197` is not substituted prospectively because the curated locus reference in the candidate pool is `FR845719.1`.

### BGC0002048 — finished truth but G3 provenance not yet cleanly resolved

- Organism: *Mycetohabitans rhizoxinica* HKI 454
- Reference: `FR687359.1`, from a complete genome assembly
- 67 candidate rows were recovered.
- The large PRJEB76713 block comes from `Inducing Novel Endosymbioses by Bacteria Implantation into Fungi` and is not treated as a clean primary-strain run without additional provenance adjudication.
- `SRR26184017` from PRJNA996741 (`BactID`) is a separate exact-name WGS/GENOMIC/RANDOM/PAIRED MiSeq row and remains potentially usable, but its biological-source provenance requires manual verification before any lock.
- Decision: hold; do not advance ahead of cleaner, higher-ranked eligible systems.

### BGC0002061 — same *Nostoc punctiforme* biological system as BGC0001479

- Organism/run universe is the same PCC 73102 system and the same seven candidate rows as BGC0001479.
- Decision: if BGC0001479 is promoted as VAL0008, BGC0002061 becomes a within-genome secondary locus and cannot count as another independent system.

### BGC0002344 — eligible-looking reserve

- Organism: *Nonomuraea coxensis* DSM 45129
- Reference: `CP068985.1`; the later publication describes the DSM 45129 genome as completely assembled into one circular chromosome.
- Candidate run: `SRR3947903`; exact strain, WGS/GENOMIC/RANDOM/PAIRED/ILLUMINA HiSeq 2000; one candidate row.
- Decision: retain as a strong independent reserve after earlier ranked clean systems; same-strain linkage between the historical Illumina run and later finished truth must be stated explicitly if promoted.

## Prospective next action

The next permitted candidate-selection action is G3 locking of `BGC0001479 / SRR30545090` as `VAL0008`, using the frozen candidate table, reference `CP001037.1`, assembly `GCA_000020025.1`, sample `SAMN43372961`, and study `PRJNA1152567`.

No reads have been downloaded for this system and no G4/G5/ranking outcome has been inspected. The selection must not be changed based on downstream results.
