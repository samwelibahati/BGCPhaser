from typer.testing import CliRunner

from bgcphaser import __version__
from bgcphaser.cli import app


runner = CliRunner()


def test_cli_version() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == f"BGCPhaser {__version__}"
