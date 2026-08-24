# Protocol amendment A1: G5 graph-walk anchor localization

Date: 2026-08-25

Status: frozen before G6, candidate enumeration, chemistry analysis, or inspection of any BGCPhaser score for a new validation system.

## Trigger

The first expanded-validation system, VAL0001, passed G1-G4. The initial G5 implementation aligned each 250-nt endpoint anchor only to individual GFA segment sequences. Neither anchor produced a qualifying whole-anchor single-segment alignment. The raw mapping evidence showed that the left anchor was split across graph segments, with two 100%-identity fragments covering 101 nt and 143 nt of the 250-nt query, and the right anchor contained a 200-nt 100%-identity fragment on another graph segment.

No BGCPhaser candidate score, truth rank, chemistry score, or reference-similarity outcome had been inspected when this implementation issue was identified.

## Methodological determination

The biological G5 criterion is unique, orientation-consistent localization of each predefined endpoint anchor in the canonical short-read assembly graph. Requiring the complete 250-nt anchor to lie inside one GFA segment adds a graph-segmentation constraint that is not part of that criterion. A valid anchor can span one or more directed GFA links while remaining uniquely localized in the graph.

The single-segment G5 implementation is therefore retained as an audit artifact but is superseded for subsequent primary validation by graph-walk localization.

## Bias control and cohort consequence

VAL0001 is designated a protocol-amendment/development sentinel and does not contribute to the primary independent-validation performance estimates after this amendment. Any graph-walk localization, candidate enumeration, or ranking result subsequently obtained for VAL0001 is descriptive/developmental only.

The amended prospective primary-validation cohort begins with VAL0002. No system is promoted into the primary cohort on the basis of a result observed under the superseded single-segment implementation.

## Frozen amended G5 rule for VAL0002 onward

The following quantities remain unchanged:

- MIBiG coordinate convention: 1-based inclusive start and 1-based exclusive end;
- left anchor: first 250 nt inside the curated locus, reference-forward;
- right anchor: last 250 nt inside the curated locus, reference-forward;
- minimap2 version: 2.31;
- minimap2 preset: `asm5`;
- minimum query coverage: 95%;
- minimum nucleotide identity: 98%.

For graph-walk localization, canonical GFA segment sequences are spelled along oriented assembler path records in `contigs.paths` and `scaffolds.paths`. Consecutive states can be spelled together only when the selected directed GFA graph contains the corresponding `L`-edge transition; the recorded GFA overlap CIGAR is applied during spelling. Identical oriented graph walks are deduplicated before mapping.

### A1.1 handling of SPAdes gapped path records

During development-sentinel execution, an assembler path contained the transition `1705+ -> 125059+`, which was absent from the directed `L`-edge graph represented by `assembly_graph_with_scaffolds.gfa`. This occurred before minimap2 anchor mapping and therefore produced no G5 biological outcome.

Assembler path records can contain scaffold jumps or other path-level adjacency that is not a sequence-spellable directed `L` edge. Such a transition is not used to infer nucleotide sequence and is not introduced into the BGCPhaser graph. Each assembler path record is instead split at every adjacent-state transition absent from the directed `L` graph. The resulting maximal contiguous directed graph-walk chunks are retained, single-state chunks included, and exact duplicate oriented chunks are deduplicated before sequence spelling and anchor mapping. No unknown bases, arbitrary gap sequence, or synthetic graph edge is inserted across a split.

This A1.1 rule is a format/topology implementation clarification and applies uniformly to all subsequent graph-walk G5 analyses, beginning before the first post-amendment primary system VAL0002.

Each 250-nt reference anchor is aligned to the spelled contiguous graph-walk sequences. Qualifying alignments may span one or more linked oriented GFA states. Equivalent alignments arising from overlapping contig/scaffold path records or reverse-complement path representations are collapsed to the same physical graph-walk localization using the reference-forward oriented state subpath and its boundary-state offsets.

Exactly one qualifying physical graph-walk localization is required for the left anchor and exactly one for the right anchor. The unique reference-forward left-anchor walk defines the left boundary state, and the unique reference-forward right-anchor walk defines the right boundary state for G6. Zero or more than one physical localization for either anchor fails G5.

Anchor length, identity threshold, coverage threshold, minimap2 preset, or graph-walk definition are not altered for an individual system after its mapping result is observed. The amended rule applies uniformly to every primary candidate from VAL0002 onward.
