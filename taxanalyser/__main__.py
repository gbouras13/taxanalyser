"""
Entrypoint for taxanlayser

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("taxanlayser.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("taxanlayser.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-o",
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="taxanlayser.out",
            show_default=True,
        ),
        click.option(
            "--config",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "-t",
            "--threads",
             help="Number of threads to use", 
             default=1, show_default=True
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="taxanlayser.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Snakemake and Snaketool pipeline to taxonomically profile ONT long read metagenomics with sourmash
    \b
    For more options, run:
    taxanlayser command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
taxanlayser run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           taxanlayser run --input [file]
Specify threads:    taxanlayser run ... --threads [threads]
Disable conda:      taxanlayser run ... --no-use-conda 
Change defaults:    taxanlayser run ... --snake-default="-k --nolock"
Add Snakemake args: taxanlayser run ... --dry-run --keep-going --touch
Specify targets:    taxanlayser run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_msg_install = """
\b
Download the host genome 
taxanlayser install ... 
\b
RUN EXAMPLES:
Database:           taxanlayser install --database [file]
"""



@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i","--input", "_input", help="Input csv with 2 columns", type=str, required=True)
@common_options
def long(_input, output, log, config, **kwargs):
    """Run taxanlayser long"""
    # Config to add or update in configfile
    merge_config = { "args": {
        "input": _input,
        "output": output,
        "config": config,
        "log": log
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_taxanlayser.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_msg_install,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ))
@click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        )
@click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
                "--conda-frontend conda"
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        )
@click.option(
            "--database",
            help="Database directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="host_db",
            show_default=True,
        )
@common_options
def install( database, log,   **kwargs):
    """Install host DB"""
    # Config to add or update in configfile
    merge_config = {
        "database": database,
        "log": log
    }
    """Install databases"""
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow','download_chm13_host.smk')),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        **kwargs)


@click.command()
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=snake_base(os.path.join("config", "config.yaml")))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(long)
# cli.add_command(paired)
cli.add_command(install)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
