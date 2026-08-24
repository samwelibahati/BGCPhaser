from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .pipeline import analyse_candidates
from .scoring import (
    CandidateScore,
    format_q12,
    rank_candidates,
    score_candidate,
)

app = typer.Typer(no_args_is_help=True)
console = Console()


def _version_callback(value: bool) -> bool:
    if value:
        typer.echo(f"BGCPhaser {__version__}")
        raise typer.Exit()
    return value


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show BGCPhaser version and exit.",
    ),
) -> None:
    """BGCPhaser command-line interface."""


def _read_candidates(path: Path) -> list[CandidateScore]:
    required = [
        "candidate_id",
        "sequence_score",
        "chemistry_c",
        "chemistry_m",
    ]
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [name for name in required if name not in fields]
        if missing:
            raise typer.BadParameter(
                "Missing input columns: " + ", ".join(missing)
            )
        results = []
        for row_number, row in enumerate(reader, start=2):
            try:
                results.append(
                    score_candidate(
                        row["candidate_id"],
                        row["sequence_score"],
                        row["chemistry_c"],
                        row["chemistry_m"],
                    )
                )
            except ValueError as exc:
                raise typer.BadParameter(
                    f"Input row {row_number}: {exc}"
                ) from exc
    if not results:
        raise typer.BadParameter("Input contains no candidates")
    return results


@app.command()
def rank(
    input_path: Path = typer.Option(
        ...,
        "--input",
        exists=True,
        readable=True,
        dir_okay=False,
        help="TSV with candidate_id, sequence_score, chemistry_c, and chemistry_m.",
    ),
    output_path: Path = typer.Option(
        ...,
        "--output",
        dir_okay=False,
        help="Output ranked TSV.",
    ),
    top_k: int = typer.Option(10, "--top-k", min=1),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite an existing output.",
    ),
) -> None:
    """Score and rank candidates using the frozen sequence-plus-chemistry model."""
    if output_path.exists() and not force:
        raise typer.BadParameter(f"Output already exists: {output_path}")

    ranked = rank_candidates(_read_candidates(input_path))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "candidate_id",
        "sequence_score",
        "chemistry_c",
        "chemistry_m",
        "chemistry_score",
        "combined_score",
        "rank",
        "status",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for item in ranked:
            writer.writerow(
                {
                    "candidate_id": item.candidate_id,
                    "sequence_score": format_q12(item.sequence_score),
                    "chemistry_c": format_q12(item.chemistry_c),
                    "chemistry_m": format_q12(item.chemistry_m),
                    "chemistry_score": format_q12(item.chemistry_score),
                    "combined_score": format_q12(item.combined_score),
                    "rank": str(item.rank) if item.rank is not None else "NA",
                    "status": item.status,
                }
            )

    table = Table(title="BGCPhaser ranked candidates")
    for column in [
        "rank",
        "candidate",
        "S_seq",
        "C",
        "M",
        "S_chem",
        "S_combined",
        "status",
    ]:
        table.add_column(column)
    for item in ranked[:top_k]:
        table.add_row(
            str(item.rank) if item.rank is not None else "NA",
            item.candidate_id,
            format_q12(item.sequence_score),
            format_q12(item.chemistry_c),
            format_q12(item.chemistry_m),
            format_q12(item.chemistry_score),
            format_q12(item.combined_score),
            item.status,
        )
    console.print(table)


@app.command()
def analyse(
    gfa: Path = typer.Option(
        ...,
        "--gfa",
        exists=True,
        readable=True,
        dir_okay=False,
        help="Assembly graph in GFA1 format.",
    ),
    candidates: Path = typer.Option(
        ...,
        "--candidates",
        exists=True,
        readable=True,
        dir_okay=False,
        help="TSV with candidate_id and oriented_walk.",
    ),
    read1: Path = typer.Option(
        ...,
        "--read1",
        exists=True,
        readable=True,
        dir_okay=False,
        help="Paired-end read 1 FASTQ.",
    ),
    read2: Path = typer.Option(
        ...,
        "--read2",
        exists=True,
        readable=True,
        dir_okay=False,
        help="Paired-end read 2 FASTQ.",
    ),
    unpaired_read: list[Path] = typer.Option(
        [],
        "--unpaired-read",
        exists=True,
        readable=True,
        dir_okay=False,
        help="Optional unpaired short-read FASTQ for R evidence. Repeat for multiple files.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        file_okay=False,
        help="New output directory. Existing directories are not overwritten.",
    ),
    insert_min: int = typer.Option(
        ...,
        "--insert-min",
        min=0,
        help="Locked empirical minimum paired-fragment length.",
    ),
    insert_max: int = typer.Option(
        ...,
        "--insert-max",
        min=0,
        help="Locked empirical maximum paired-fragment length; also defines L transition-context flank.",
    ),
    antismash_databases: Path = typer.Option(
        ...,
        "--antismash-databases",
        exists=True,
        readable=True,
        file_okay=False,
        help="antiSMASH database directory.",
    ),
    contigs_path: Optional[Path] = typer.Option(
        None,
        "--contigs-path",
        exists=True,
        readable=True,
        dir_okay=False,
        help="Optional SPAdes contigs.paths for H.",
    ),
    scaffolds_path: Optional[Path] = typer.Option(
        None,
        "--scaffolds-path",
        exists=True,
        readable=True,
        dir_okay=False,
        help="Optional SPAdes scaffolds.paths for H.",
    ),
    depth_tag: Optional[str] = typer.Option(
        None,
        "--depth-tag",
        help="GFA depth tag for V, for example DP or dp. Omit to make V NA.",
    ),
    minimap2_executable: str = typer.Option(
        "minimap2",
        "--minimap2",
        help="minimap2 executable.",
    ),
    antismash_executable: str = typer.Option(
        "antismash",
        "--antismash",
        help="antiSMASH executable.",
    ),
    minimap2_threads: int = typer.Option(1, "--minimap2-threads", min=1),
    antismash_cpus: int = typer.Option(1, "--antismash-cpus", min=1),
) -> None:
    """Run intrinsic sequencing and chemistry evidence and rank candidate BGC reconstructions."""
    assembler_paths = [
        path for path in (contigs_path, scaffolds_path) if path is not None
    ]
    try:
        result = analyse_candidates(
            gfa_path=gfa,
            candidate_path=candidates,
            read1_path=read1,
            read2_path=read2,
            output_dir=output_dir,
            insert_min=insert_min,
            insert_max=insert_max,
            antismash_database_root=antismash_databases,
            unpaired_read_paths=unpaired_read,
            assembler_paths=assembler_paths,
            depth_tag=depth_tag,
            minimap2_executable=minimap2_executable,
            antismash_executable=antismash_executable,
            minimap2_threads=minimap2_threads,
            antismash_cpus=antismash_cpus,
        )
    except (
        FileNotFoundError,
        FileExistsError,
        PermissionError,
        RuntimeError,
        ValueError,
    ) as exc:
        console.print(f"[red]BGCPhaser analyse failed:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    console.print("[green]BGCPhaser analyse complete[/green]")
    console.print(f"Candidates: {result.candidate_count}")
    console.print(f"Ranked: {result.ranked_count}")
    console.print(f"Results: {result.ranked_path}")


if __name__ == "__main__":
    app()
