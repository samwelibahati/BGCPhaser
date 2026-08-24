# BGCPhaser

**BGCPhaser** ranks alternative modular biosynthetic gene cluster (BGC) reconstructions from nucleotide assembly graphs by integrating intrinsic sequencing evidence with biosynthetic-architecture evidence.

Version: **0.1.0**

## Scope

BGCPhaser is designed for ambiguous modular bacterial BGC reconstructions, particularly nonribosomal peptide synthetase (NRPS), type I polyketide synthase (T1PKS), and hybrid NRPS-T1PKS loci represented by multiple assembler-bounded graph walks.

The central design principle is analytical separation: reference BGC sequences, MIBiG similarity and known truth sequences are not used to compute candidate scores or ranks. Reference information may be used separately to define benchmark loci and to evaluate a frozen ranking after scoring.

## Evidence model

For each candidate, BGCPhaser can integrate four sequencing-evidence components:

- `R`: direct read support across candidate graph transitions;
- `L`: paired-fragment physical linkage across transitions;
- `V`: copy-aware graph-segment coverage consistency;
- `1-H`: concordance with assembler-native path adjacencies.

It also derives two biosynthetic-architecture components from antiSMASH modular NRPS/PKS annotations:

- `C`: intrinsic loader/substrate support among assessable modules;
- `M`: module completeness.

The frozen chemistry score is

```text
S_chem = q12(0.5*C + 0.5*M)
```

and the frozen combined score is

```text
S_combined = q12(0.5*S_seq + 0.5*S_chem)
```

where `q12` denotes 12-decimal `ROUND_HALF_EVEN` quantization. Combined ranking is complete-case, and exact q12 score ties use standard competition ranking.

## Requirements

- Python 3.11 or later
- minimap2
- antiSMASH with its databases for chemistry-enabled analysis

The benchmarked workflow used minimap2 2.31-r1302 and antiSMASH 8.0.4. Other external-tool versions have not yet been established as benchmark-equivalent.

## Installation

```bash
git clone https://github.com/samwelibahati/BGCPhaser.git
cd BGCPhaser
python -m pip install -e .
```

Check the command-line interface:

```bash
bgcphaser --version
bgcphaser --help
```

## Candidate input

The `analyse` workflow accepts a tab-separated candidate table containing at minimum:

```text
candidate_id    oriented_walk
```

An oriented walk contains graph-segment identifiers with explicit orientations, for example:

```text
candidate_001    101+->205-->307+
```

Candidate-locus localization and benchmark truth evaluation are intentionally separate from intrinsic ranking.

## Main commands

Run end-to-end evidence extraction and ranking:

```bash
bgcphaser analyse --help
```

Rank candidates from already-computed features:

```bash
bgcphaser rank --help
```

## Validation status

BGCPhaser v0.1.0 was frozen against a real-data benchmark workflow before public release. The BGRR0074 regression fixture reproduced the locked sequence, chemistry and combined scores and competition rank. The principal ambiguous benchmark, BGRR0055, was used to evaluate prioritization among a large assembler-bounded candidate population. Independent prospective systems were screened under predeclared eligibility criteria without changing the frozen method.

These benchmarks define the empirical scope of v0.1.0; they are not evidence that all BGC classes, assemblers or sequencing designs have equivalent performance.

## Tests

Run the package tests with:

```bash
python -m pytest -q
```

## Reproducibility

The repository includes source code, automated tests, environment metadata and small test fixtures. Public sequencing reads and curated BGC records used in the manuscript remain available from their source repositories (NCBI SRA and MIBiG) and are not redistributed here.

## Citation

Citation metadata are provided in `CITATION.cff`. Please cite the archived software release DOI once available, together with the associated BGCPhaser manuscript when published.

## License

BGCPhaser is released under the MIT License. See `LICENSE`.
