#!/usr/bin/python3
# -*-coding:Latin-1 -*


"""
Licence header
--------------
dchimer : python3 main program of d-chimer pipeline.
Copyright (C) 2022, INSTITUT PASTEUR DE LA GUYANE
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.
date: 19.04.2022
--
Author : Sourakhata Tirera <stirera@pasteur-cayenne.fr>
Author : Jean-Marc Frigerio <jean-marc.frigerio@inrae.fr>
Author : Alix de Thoisy <alixdet@protonmail.com>
"""

import os
import sys

from .dchimer_methods import call_blastx_and_filter, call_blastn_and_filter
from .dchimer_logging import setup_logging, get_logger
from .dchimer_config_validator import (
    load_and_validate_config,
    validate_database_paths,
    validate_blast_executable
)
from .dchimer_input_validator import validate_fasta_file

import argparse


def opts_and_args():
    """d-chimer arguments and options"""
    description = 'run a blast program and filter its outputs'
    usage = 'usage : dchimer [-h] -p blastProgram -L True -f fastafile [-m] max_loops'
    parser = argparse.ArgumentParser(description=description, epilog=usage)

    parser.add_argument('-p',
                        '--program',
                        type=str,
                        help='blastProgram (required) : blastn or blastx',
                        dest='blastprogram',
                        metavar='blastprogram',
                        required=True,
                        default='blastn')

    parser.add_argument('-L',
                        '--local',
                        type=bool,
                        help="""use local machine BLAST+ programs or Biopython
                                embedded ones (required) :
                                True = use local AND FALSE = use Biopython BLAST+ (default : True)""",
                        dest='local',
                        metavar='Tue | False',
                        required=True,
                        default=True)

    parser.add_argument('-f',
                        '--fasta',
                        type=str,
                        help='input fasta_file (required)',
                        dest='fastafile',
                        metavar='fastafile',
                        required=True,
                        default='')

    parser.add_argument('-m',
                        '--max_recursive_loops',
                        type=int,
                        help="""maximum number of times to to
                                process uncovered zones fasta""",
                        dest='max_loops',
                        metavar='max_loops',
                        required=False,
                        default=10000)
    args = parser.parse_args()
    return args


def main():
    """Main entry point for d-chimer pipeline.
    
    Validates configuration and input, then runs BLAST and filtering
    until max_loops cycles reached or process completes.
    """
    # Initialize logging
    logger = setup_logging()
    logger.info("Starting d-chimer pipeline")

    try:
        # Parse arguments
        args = opts_and_args()

        program = args.blastprogram
        local = args.local
        qfile = args.fastafile
        cpt = 0
        cpt_max = args.max_loops

        # Validate program type
        if program not in ['blastn', 'blastx']:
            logger.error(f"Invalid BLAST program: {program}")
            logger.error("Allowed values: blastn or blastx")
            sys.exit(1)

        logger.info(f"BLAST program: {program}")
        logger.info(f"Local BLAST: {local}")
        logger.info(f"Input file: {qfile}")
        logger.info(f"Maximum cycles: {cpt_max}")

        # Load and validate configuration
        logger.info("Validating configuration...")
        config_path = os.path.join(
            os.path.dirname(__file__),
            "dchimer_config.yaml"
        )
        config = load_and_validate_config(config_path)
        logger.info("Configuration loaded successfully")

        # Validate database paths
        logger.info("Validating database paths...")
        validate_database_paths(config)
        logger.info("Database paths validated")

        # Validate BLAST executables if using local BLAST
        if local:
            logger.info("Validating BLAST installation...")
            validate_blast_executable(config['blast_path'])
            logger.info("BLAST installation validated")

        # Validate input FASTA file
        logger.info(f"Validating input FASTA file: {qfile}")
        num_seqs, file_size = validate_fasta_file(qfile)
        logger.info(f"Input file valid: {num_seqs} sequences, {file_size:.2f} KB")

        # Check output directory doesn't already exist
        root = qfile.split(".")
        if program == 'blastn':
            output_dir = root[0] + "_bn_out"
        else:
            output_dir = root[0] + "_bx_out"

        if os.path.isdir(output_dir):
            logger.error(f"Output directory already exists: {output_dir}")
            logger.error("Please remove it before running d-chimer")
            sys.exit(1)

        logger.info(f"Output directory: {output_dir}")
        logger.info("All validation passed. Starting BLAST and filtering...")

        # Run appropriate BLAST pipeline
        if program == 'blastn':
            call_blastn_and_filter(program, local, qfile, cpt, cpt_max)
        else:  # blastx
            call_blastx_and_filter(program, local, qfile, cpt, cpt_max)

        logger.info("d-chimer pipeline completed successfully")

    except (FileNotFoundError, PermissionError, ValueError) as e:
        logger.error(f"Validation error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Pipeline interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
