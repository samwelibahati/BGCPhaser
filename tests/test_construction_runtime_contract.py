from decimal import Decimal

from bgcphaser.construction import (
    candidate_transition_support,
    construct_transition_classes,
    spell_candidate,
)
from bgcphaser.gfa import read_gfa
from bgcphaser.tools import (
    FORBIDDEN_ANTISMASH_OPTIONS,
    antismash_command,
    minimap2_command,
)


def test_validated_construction_and_runtime_contract(tmp_path):
    gfa = tmp_path / "toy.gfa"
    gfa.write_text(
        "S\tA\t" + ("A" * 200) + "\n"
        + "S\tB\t" + ("C" * 200) + "\n"
        + "L\tA\t+\tB\t+\t127M\n",
        encoding="ascii",
    )
    graph = read_gfa(gfa)
    path = [("A", "+"), ("B", "+")]
    spelling = spell_candidate(graph, path)
    assert len(spelling.sequence) == 273
    assert spelling.junctions == (200,)
    assert spelling.overlaps == (127,)

    classes = construct_transition_classes(
        graph,
        {"candidate_1": path},
        target_flank=160,
        minimum_side=30,
        class_prefix="TOY",
    )
    assert len(classes.r_sequences) == 1
    assert len(classes.l_sequences) == 1
    r_class_id = next(iter(classes.r_sequences))
    l_class_id = next(iter(classes.l_sequences))
    assert len(classes.r_sequences[r_class_id]) == 187
    assert len(classes.l_sequences[l_class_id]) == 233
    assert classes.l_junctions[l_class_id] == 160
    assert len(classes.occurrences) == 1
    occurrence = classes.occurrences[0]
    assert occurrence.r_class_id == r_class_id
    assert occurrence.l_class_id == l_class_id
    assert candidate_transition_support(
        classes.occurrences,
        candidate_id="candidate_1",
        evidence="R",
        class_support={r_class_id: 5},
    ) == Decimal("1")
    assert candidate_transition_support(
        classes.occurrences,
        candidate_id="candidate_1",
        evidence="L",
        class_support={l_class_id: 0},
    ) == Decimal("0")

    mm2 = minimap2_command(
        executable="minimap2",
        target_fasta="classes.fa",
        reads=["reads_R1.fastq.gz", "reads_R2.fastq.gz"],
        threads=8,
    )
    assert mm2 == [
        "minimap2", "-a", "-x", "sr", "--secondary=yes", "-N", "1000",
        "-p", "0", "-t", "8", "classes.fa", "reads_R1.fastq.gz",
        "reads_R2.fastq.gz",
    ]

    antismash = antismash_command(
        executable="antismash",
        input_fasta="candidate.fa",
        output_dir="candidate.out",
        output_basename="candidate_1",
        database_root="databases",
        logfile="candidate.log",
        cpus=4,
    )
    assert "--minimal" in antismash
    assert "--enable-nrps-pks" in antismash
    assert not (FORBIDDEN_ANTISMASH_OPTIONS & set(antismash))
    assert "--cc-mibig" not in antismash
    assert "--reuse-results" not in antismash
