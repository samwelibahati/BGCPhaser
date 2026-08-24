from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
from typing import Sequence


FORBIDDEN_ANTISMASH_OPTIONS = {
    "--cb-general",
    "--cb-knownclusters",
    "--cb-subclusters",
    "--cc-mibig",
    "--cc-custom-dbs",
    "--reuse-results",
    "--sideload",
    "--sideload-by-cds",
    "--sideload-simple",
}


def resolve_executable(
    executable: str | Path,
) -> Path:
    """
    Resolve an explicit executable path or a command on PATH.
    """
    text = str(
        executable
    )

    expanded = Path(
        text
    ).expanduser()

    if (
        expanded.is_absolute()
        or "/" in text
    ):
        resolved = expanded

    else:
        found = shutil.which(
            text
        )

        if found is None:
            raise FileNotFoundError(
                f"Executable not found: "
                f"{text}"
            )

        resolved = Path(
            found
        )

    if not resolved.is_file():
        raise FileNotFoundError(
            f"Executable not found: "
            f"{resolved}"
        )

    if not os.access(
        resolved,
        os.X_OK,
    ):
        raise PermissionError(
            f"Not executable: "
            f"{resolved}"
        )

    return resolved.resolve()


def minimap2_command(
    *,
    executable: str | Path,
    target_fasta: str | Path,
    reads: Sequence[
        str | Path
    ],
    threads: int,
) -> list[str]:
    """
    Build the validated primitive-evidence minimap2 command.

    Primary and secondary non-supplementary alignments
    remain available for strict AS competition.
    """
    if threads < 1:
        raise ValueError(
            "threads must be >= 1"
        )

    if not reads:
        raise ValueError(
            "At least one read file "
            "is required"
        )

    return [
        str(executable),
        "-a",
        "-x",
        "sr",
        "--secondary=yes",
        "-N",
        "1000",
        "-p",
        "0",
        "-t",
        str(threads),
        str(target_fasta),
        *[
            str(path)
            for path in reads
        ],
    ]


def run_minimap2(
    *,
    executable: str | Path,
    target_fasta: str | Path,
    reads: Sequence[
        str | Path
    ],
    sam_path: str | Path,
    log_path: str | Path,
    threads: int = 1,
) -> Path:
    """
    Execute minimap2 without shell interpolation and
    atomically publish the completed SAM file.
    """
    executable_path = (
        resolve_executable(
            executable
        )
    )

    target = Path(
        target_fasta
    )

    if not target.is_file():
        raise FileNotFoundError(
            f"Target FASTA absent: "
            f"{target}"
        )

    read_paths = [
        Path(path)
        for path in reads
    ]

    for read_path in read_paths:
        if not read_path.is_file():
            raise FileNotFoundError(
                f"Read file absent: "
                f"{read_path}"
            )

    output = Path(
        sam_path
    )

    log = Path(
        log_path
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    log.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = Path(
        str(output)
        + ".tmp"
    )

    if temporary.exists():
        raise FileExistsError(
            f"Temporary SAM exists: "
            f"{temporary}"
        )

    command = minimap2_command(
        executable=(
            executable_path
        ),
        target_fasta=target,
        reads=read_paths,
        threads=threads,
    )

    try:
        with (
            temporary.open(
                "w",
                encoding="utf-8",
            ) as stdout_handle,
            log.open(
                "w",
                encoding="utf-8",
            ) as stderr_handle,
        ):
            completed = subprocess.run(
                command,
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
                check=False,
            )

        if completed.returncode != 0:
            temporary.unlink(
                missing_ok=True
            )

            raise RuntimeError(
                "minimap2 failed with "
                f"exit code "
                f"{completed.returncode}; "
                f"see {log}"
            )

        os.replace(
            temporary,
            output,
        )

    finally:
        if temporary.exists():
            temporary.unlink()

    return output


def antismash_command(
    *,
    executable: str | Path,
    input_fasta: str | Path,
    output_dir: str | Path,
    output_basename: str,
    database_root: str | Path,
    logfile: str | Path,
    cpus: int,
) -> list[str]:
    """
    Build the validated intrinsic-chemistry antiSMASH
    command.

    No ClusterBlast, ClusterCompare, MIBiG comparison,
    result reuse, or sideloading option is admitted.
    """
    if cpus < 1:
        raise ValueError(
            "cpus must be >= 1"
        )

    if not output_basename:
        raise ValueError(
            "output_basename is required"
        )

    command = [
        str(executable),
        str(input_fasta),
        "--taxon",
        "bacteria",
        "--genefinding-tool",
        "prodigal",
        "--databases",
        str(database_root),
        "--minimal",
        "--enable-nrps-pks",
        "--no-enable-html",
        "--no-zip-output",
        "--no-region-gbks",
        "--no-summary-gbk",
        "--abort-on-invalid-records",
        "--cpus",
        str(cpus),
        "--output-dir",
        str(output_dir),
        "--output-basename",
        output_basename,
        "--logfile",
        str(logfile),
    ]

    if (
        FORBIDDEN_ANTISMASH_OPTIONS
        & set(command)
    ):
        raise RuntimeError(
            "Reference-comparison or "
            "prior-result option entered "
            "the chemistry command"
        )

    return command


def run_antismash(
    *,
    executable: str | Path,
    input_fasta: str | Path,
    output_dir: str | Path,
    output_basename: str,
    database_root: str | Path,
    cpus: int = 1,
) -> Path:
    """
    Execute intrinsic NRPS/PKS antiSMASH annotation.

    The final directory is published only after a successful
    execution with the expected JSON result.
    """
    executable_path = (
        resolve_executable(
            executable
        )
    )

    input_path = Path(
        input_fasta
    )

    database_path = Path(
        database_root
    )

    final_dir = Path(
        output_dir
    )

    staging_dir = Path(
        str(final_dir)
        + ".staging"
    )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input FASTA absent: "
            f"{input_path}"
        )

    if not database_path.is_dir():
        raise FileNotFoundError(
            f"antiSMASH database "
            f"directory absent: "
            f"{database_path}"
        )

    if final_dir.exists():
        raise FileExistsError(
            f"Output directory exists: "
            f"{final_dir}"
        )

    if staging_dir.exists():
        raise FileExistsError(
            f"Staging directory exists: "
            f"{staging_dir}"
        )

    staging_dir.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    staging_dir.mkdir()

    logfile = (
        staging_dir
        / "antismash.console.log"
    )

    command = antismash_command(
        executable=(
            executable_path
        ),
        input_fasta=input_path,
        output_dir=staging_dir,
        output_basename=(
            output_basename
        ),
        database_root=(
            database_path
        ),
        logfile=logfile,
        cpus=cpus,
    )

    runtime_env = os.environ.copy()

    executable_bin = str(
        executable_path.parent
    )

    inherited_path = (
        runtime_env.get(
            "PATH",
            "",
        )
    )

    runtime_env["PATH"] = (
        executable_bin
        + (
            os.pathsep
            + inherited_path
            if inherited_path
            else ""
        )
    )

    completed = subprocess.run(
        command,
        check=False,
        env=runtime_env,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "antiSMASH failed with "
            f"exit code "
            f"{completed.returncode}; "
            f"staging retained at "
            f"{staging_dir}"
        )

    expected_json = (
        staging_dir
        / f"{output_basename}.json"
    )

    if not expected_json.is_file():
        raise RuntimeError(
            "antiSMASH completed without "
            f"expected JSON: "
            f"{expected_json}"
        )

    os.replace(
        staging_dir,
        final_dir,
    )

    return (
        final_dir
        / f"{output_basename}.json"
    )
