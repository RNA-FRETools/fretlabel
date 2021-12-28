#!/usr/bin/env python3

import subprocess
import os
import argparse

from fretlabel import metadata

package_directory = os.path.dirname(os.path.abspath(__file__))


def solvate():
    subprocess.call(os.path.join(package_directory, "skripts", "solvate.sh"))


def single_run():
    subprocess.call(os.path.join(package_directory, "skripts", "single_run.sh"))


def continue_run():
    subprocess.call(os.path.join(package_directory, "skripts", "continue_run.sh"))


def multi_run():
    subprocess.call(os.path.join(package_directory, "skripts", "multi_run.sh"))


def resp_fit():
    subprocess.call(os.path.join(package_directory, "skripts", "resp_fit.sh"))


def fretlabel():
    parser = argparse.ArgumentParser(
        description="Interactively label nucleic acids with fluorophores"
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + str(metadata["Version"])
    )
    parser.add_argument(
        "--path",
        action="version",
        version=f"package directory: {package_directory}",
        help="Show package directory",
    )
    parser.parse_args()
