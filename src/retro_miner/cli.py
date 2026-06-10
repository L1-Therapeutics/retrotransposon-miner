from __future__ import annotations

import click


@click.group()
def cli() -> None:
    """Retrotransposon miner command-line interface."""


@cli.command("check-env")
def check_env() -> None:
    """Print basic environment status."""
    click.echo("retrotransposon-miner CLI is installed.")
    click.echo("Run scripts/validate_environment.sh for full validation.")


if __name__ == "__main__":
    cli()
