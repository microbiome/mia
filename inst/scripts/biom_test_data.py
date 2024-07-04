#!/usr/bin/env python3
# Copyright (c) 2024, Moritz E. Beber.


# /// script
# requires-python = ">=3.8"
# dependencies = [
#   "numpy ~=1.25",
#   "biom-format ~=2.1",
# ]
# ///


"""Generate BIOM files for testing purposes."""


import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from biom import Table
from biom.util import biom_open


logger = logging.getLogger()


def generate_complete_table() -> Table:
    """Generate a complete BIOM table with full taxonomy information."""
    return Table(
        np.asarray(
            [
                [100, 120],
                [90, 110],
                [80, 100],
                [70, 90],
                [40, 50],
                [30, 40],
            ]
        ),
        observation_ids=[
            "NCBI:txid2",
            "NCBI:txid1783272",
            "NCBI:txid186801",
            "NCBI:txid3085636",
            "NCBI:txid186803",
            "NCBI:txid2603322",
        ],
        sample_ids=["A", "B"],
        observation_metadata=[
            {"taxonomy": ["Bacteria"]},
            {"taxonomy": ["Bacteria", "Terrabacteria group"]},
            {"taxonomy": ["Bacteria", "Terrabacteria group", "Clostridia"]},
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                ]
            },
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                    "Lachnospiraceae",
                ]
            },
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                    "Vallitaleaceae",
                ]
            },
        ],
        # FIXME (Moritz): This is not working as expected, see https://github.com/biocore/biom-format/issues/971.
        # observation_group_metadata={"ranks": ("list", ["Superkingdom", "Clade", "Class", "Order", "Family"])},
        observation_group_metadata={
            "ranks": ("csv", "Superkingdom;Clade;Class;Order;Family")
        },
        type="Taxon table",
    )


def generate_simple_table() -> Table:
    """Generate a simple BIOM table without any taxonomy information."""
    return Table(
        np.asarray(
            [
                [100, 120],
                [90, 110],
                [80, 100],
                [70, 90],
                [40, 50],
                [30, 40],
            ]
        ),
        observation_ids=[
            "NCBI:txid2",
            "NCBI:txid1783272",
            "NCBI:txid186801",
            "NCBI:txid3085636",
            "NCBI:txid186803",
            "NCBI:txid2603322",
        ],
        sample_ids=["A", "B"],
        type="Taxon table",
    )


def generate_less_ranks_table() -> Table:
    """Generate a BIOM table with less ranks than levels in the lineage."""
    return Table(
        np.asarray(
            [
                [100, 120],
                [90, 110],
                [80, 100],
                [70, 90],
                [40, 50],
                [30, 40],
            ]
        ),
        observation_ids=[
            "NCBI:txid2",
            "NCBI:txid1783272",
            "NCBI:txid186801",
            "NCBI:txid3085636",
            "NCBI:txid186803",
            "NCBI:txid2603322",
        ],
        sample_ids=["A", "B"],
        observation_metadata=[
            {"taxonomy": ["Bacteria"]},
            {"taxonomy": ["Bacteria", "Terrabacteria group"]},
            {"taxonomy": ["Bacteria", "Terrabacteria group", "Clostridia"]},
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                ]
            },
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                    "Lachnospiraceae",
                ]
            },
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                    "Vallitaleaceae",
                ]
            },
        ],
        observation_group_metadata={"ranks": ("csv", "Superkingdom;Class;Family")},
        type="Taxon table",
    )


def generate_more_ranks_table() -> Table:
    """Generate a BIOM table with more ranks than levels in the lineage."""
    return Table(
        np.asarray(
            [
                [100, 120],
                [90, 110],
                [80, 100],
                [70, 90],
                [40, 50],
                [30, 40],
            ]
        ),
        observation_ids=[
            "NCBI:txid2",
            "NCBI:txid1783272",
            "NCBI:txid186801",
            "NCBI:txid3085636",
            "NCBI:txid186803",
            "NCBI:txid2603322",
        ],
        sample_ids=["A", "B"],
        observation_metadata=[
            {"taxonomy": ["Bacteria"]},
            {"taxonomy": ["Bacteria", "Terrabacteria group"]},
            {"taxonomy": ["Bacteria", "Terrabacteria group", "Clostridia"]},
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                ]
            },
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                    "Lachnospiraceae",
                ]
            },
            {
                "taxonomy": [
                    "Bacteria",
                    "Terrabacteria group",
                    "Clostridia",
                    "Lachnospirales",
                    "Vallitaleaceae",
                ]
            },
        ],
        observation_group_metadata={
            "ranks": ("csv", "Superkingdom;Clade;Class;Order;Family;Genus;Species")
        },
        type="Taxon table",
    )


def write_biom(table: Table, path: Path) -> None:
    """Write a BIOM table to a file."""
    with biom_open(path, permission="w") as handle:
        table.to_hdf5(
            handle,
            generated_by="taxpasta==0.7.0",
            creation_date=datetime.now(tz=timezone.utc),
        )


def main(data_dir: Path) -> None:
    """Execute the generation of test BIOM files."""
    logger.info("Generate a complete BIOM table with full taxonomy information.")
    write_biom(generate_complete_table(), data_dir / "complete.biom")
    logger.info("Generate a simple BIOM table without any taxonomy information.")
    write_biom(generate_simple_table(), data_dir / "simple.biom")
    logger.info("Generate a BIOM table with less ranks than levels in the lineage.")
    write_biom(generate_less_ranks_table(), data_dir / "less_ranks.biom")
    logger.info("Generate a BIOM table with more ranks than levels in the lineage.")
    write_biom(generate_more_ranks_table(), data_dir / "more_ranks.biom")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    main(Path(__file__).parents[1] / "extdata" / "testdata")

