# BGCPhaser expanded validation protocol v1

Status: frozen before scoring any new validation benchmark

## Objective

Evaluate whether the frozen BGCPhaser v0.1.0 intrinsic scoring model generalizes across independently selected ambiguous modular biosynthetic gene cluster (BGC) assembly graphs.

The primary validation question is whether the independently known sequence-level truth is enriched toward the top of the BGCPhaser ranking across multiple real biological systems.

## Frozen software

- BGCPhaser version: 0.1.0
- GitHub tag: v0.1.0
- Zenodo DOI: 10.5281/zenodo.22079283
- Scoring weights, missing-evidence rules, q12 arithmetic and competition ranking are fixed.
- New benchmark outcomes cannot be used to alter the v0.1.0 score before primary validation is complete.

## Source universe

Candidate systems are drawn from public records with all of the following:

1. a curated modular NRPS, type I PKS, or hybrid NRPS-PKS BGC;
2. a complete or high-confidence finished bacterial reference sequence that contains the curated locus;
3. traceable real paired-end Illumina whole-genome sequencing reads from the same isolate or a defensibly matched biological source;
4. raw read files retrievable from SRA/ENA;
5. sufficient metadata to establish provenance before short-read assembly.

MIBiG v4.0 is the primary curated BGC source. SRA/ENA provide raw sequencing data and sequencing-layout metadata.

## Exclusion of development systems from primary external validation

The previously analysed systems BGRR0055, BGRR0074, BGRR0091, BGRR0079, BGRR0057, BGRR0081 and BGRR0076 are retained as development, control, reproducibility or prior prospective-screen systems. They do not contribute to the primary independent-validation performance estimates.

## Prospective eligibility gates

Eligibility is evaluated in the following order. Screening stops at the first failed gate.

### G1. Curated modular BGC

The locus is classified as NRPS, type I PKS, or hybrid NRPS-PKS and has a traceable curated BGC record.

### G2. Finished sequence-level truth

A complete or high-confidence finished reference sequence contains the BGC locus and supports an unambiguous truth interval.

### G3. Real paired-end Illumina WGS

A traceable paired-end Illumina whole-genome run exists for the same isolate or a predeclared defensibly matched biological source. Single-end-only, mate-incomplete, amplicon, metagenomic and targeted-enrichment datasets fail this gate.

### G4. Reproducible short-read assembly

The locked read-processing and SPAdes workflow produces the canonical GFA and assembler-path records required for candidate reconstruction and intrinsic evidence extraction.

### G5. Anchor localization

Both predefined BGC endpoint anchors localize uniquely and orientation-consistently in the canonical short-read assembly graph.

### G6. Directed graph connectivity

A directed anchor-to-anchor connection exists in the oriented assembly graph.

### G7. Ambiguous assembler-bounded candidate population

Assembler-bounded enumeration terminates below the predefined search ceiling and returns at least two distinct complete candidate reconstructions. Systems with one candidate are retained as nonambiguous controls but do not enter the primary ranking-performance cohort.

### G8. Truth represented in candidate population

At least one enumerated candidate satisfies the predeclared full-window truth criterion. If the true reconstruction is absent from the candidate population, the system is unsuitable for evaluating candidate ranking and is reported separately as a candidate-generation failure.

### G9. Intrinsic scoreability

The candidate population supports the frozen BGCPhaser v0.1.0 sequence and chemistry calculations required for combined ranking.

## Outcome masking

Eligibility through G9 is determined without using BGCPhaser truth rank, MIBiG similarity rank or any outcome derived from the frozen candidate scores. Candidate scores and ranks are frozen before truth identity is joined.

## Primary outcomes

For every independently eligible ambiguous benchmark:

- candidate count;
- truth competition rank;
- truth exact-score tie size;
- normalized truth rank = truth rank / candidate count;
- reciprocal rank = 1 / truth rank;
- top-1 recovery;
- top-5 recovery;
- top-10 recovery;
- top-1%, top-5% and top-10% recovery.

Across the independent benchmark cohort:

- top-k recovery proportions;
- median truth rank;
- median normalized truth rank;
- mean reciprocal rank;
- bootstrap confidence intervals;
- performance stratified by BGC class and candidate-space size where sample size permits.

## Predeclared comparative analyses

The same eligible candidate populations are evaluated with:

1. frozen BGCPhaser v0.1.0 combined score;
2. sequence-evidence composite alone;
3. chemistry composite alone;
4. assembler-path concordance baseline;
5. coverage-consistency baseline;
6. direct read-junction support baseline;
7. paired-fragment linkage baseline;
8. post-ranking MIBiG/reference-similarity ranking;
9. random-ranking expectation;
10. biosyntheticSPAdes-derived reconstruction comparison where the output is technically commensurable.

## Predeclared ablation analyses

The frozen component values are used to evaluate leave-one-component-out ranking variants for R, L, V, H, C and M. These analyses evaluate contribution of evidence components; they do not redefine the primary v0.1.0 model.

## Computational performance

For eligible systems, record:

- candidate count;
- wall-clock time;
- CPU time where available;
- peak resident memory;
- time spent in graph/candidate processing;
- time spent in sequencing-evidence extraction;
- time spent in antiSMASH chemistry analysis;
- total output size.

## Reporting principle

All screened systems remain in the registry, including failures. Failure stage and reason are reported explicitly. No system is removed because BGCPhaser assigns the truth a low rank.
