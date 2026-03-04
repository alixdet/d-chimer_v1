#!/usr/bin/python3
# -*-coding:utf-8 -*

"""
Configuration validation for d-chimer.

Validates that required configuration parameters exist and are accessible.
"""

import os
import sys
from pathlib import Path
import yaml

def load_and_validate_config(config_path):
    """
    Load and validate d-chimer configuration.

    Args:
        config_path: path to dchimer_config.yaml

    Returns:
        config: validated configuration dictionary

    Raises:
        FileNotFoundError: if config file or referenced paths don't exist
        ValueError: if configuration is invalid
    """
    # Check config file exists
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    # Load config
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML in configuration file: {e}")
    except Exception as e:
        raise ValueError(f"Error reading configuration file: {e}")

    # Validate required sections
    required_sections = {
        'blastn_parameters': ['dbpath_nt', 'nb_threads_bn', 'evalue_nt'],
        'blastx_parameters': ['dbpath_vrl', 'dbpath_nr', 'nb_threads_bx', 
                             'evalue_vir', 'evalue_nr'],
        'filter_blastn_parameters': ['d', 'l'],
        'filter_blastx_parameters': ['d', 'l'],
        'add_taxo_parameters': ['tax_lineages_file'],
    }

    for section, params in required_sections.items():
        if section not in config:
            raise ValueError(f"Missing configuration section: '{section}'")
        
        for param in params:
            if param not in config[section]:
                raise ValueError(
                    f"Missing parameter '{param}' in section '{section}'"
                )

    # Validate blast_path
    if 'blast_path' not in config:
        raise ValueError("Missing configuration: 'blast_path'")

    return config


def validate_database_paths(config):
    """
    Validate that BLAST database paths exist and are accessible.

    Args:
        config: configuration dictionary from load_and_validate_config()

    Raises:
        FileNotFoundError: if database files don't exist
    """
    db_paths = [
        ('BLASTn NT database', config['blastn_parameters']['dbpath_nt']),
        ('BLASTx viral database', config['blastx_parameters']['dbpath_vrl']),
        ('BLASTx NR database', config['blastx_parameters']['dbpath_nr']),
        ('Taxonomy lineage file', config['add_taxo_parameters']['tax_lineages_file']),
    ]

    missing = []
    for db_name, db_path in db_paths:
        # For BLAST databases, check if any of the expected index files exist
        if 'dbpath_' in str(db_paths):
            # Check main file or with common extensions
            if not (os.path.isfile(db_path) or
                    os.path.isfile(f"{db_path}.nin") or
                    os.path.isfile(f"{db_path}.pin")):
                missing.append(f"{db_name}: {db_path}")
        else:
            # For regular files, just check existence
            if not os.path.isfile(db_path):
                missing.append(f"{db_name}: {db_path}")

    if missing:
        raise FileNotFoundError(
            f"Required database files not found:\n" +
            "\n".join(f"  - {m}" for m in missing)
        )


def validate_blast_executable(blast_path):
    """
    Validate that BLAST executables are in the specified path.

    Args:
        blast_path: path to BLAST bin directory

    Raises:
        FileNotFoundError: if BLAST executables not found
    """
    required_executables = ['blastn', 'blastx']
    missing = []

    for exe in required_executables:
        exe_path = os.path.join(blast_path, exe)
        if not os.path.isfile(exe_path) or not os.access(exe_path, os.X_OK):
            missing.append(f"{exe} in {blast_path}")

    if missing:
        raise FileNotFoundError(
            f"BLAST executables not found or not executable:\n" +
            "\n".join(f"  - {m}" for m in missing)
        )
