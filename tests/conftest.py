"""Pytest configuration and fixtures for d-chimer tests."""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_fasta_file(temp_dir):
    """Create a sample FASTA file for testing."""
    fasta_content = """>seq1
ATGCATGCATGCATGC
>seq2
GCTAGCTAGCTAGCTA
>seq3
AT"""
    
    fasta_path = os.path.join(temp_dir, "sample.fasta")
    with open(fasta_path, "w") as f:
        f.write(fasta_content)
    
    return fasta_path


@pytest.fixture
def sample_csv_blastn(temp_dir):
    """Create a sample BLAST output CSV file for BLASTn testing."""
    csv_content = """query1	gi|123456	95.5	100	4	1	1	100	1	100	1.5e-50	185	562	Archaea	16	ATGCATGCATGCATGC	ATGCATGCATGCATGC
query1	gi|789012	92.0	100	8	0	1	100	5	104	3.2e-40	155	562	Bacteria	16	ATGCATGCATGCATGC	ATGAATGCATGCATGC
query2	gi|345678	98.5	50	1	0	51	100	201	250	2.1e-25	98	562	Archaea	100	GCTAGCTAGCTAGCTA	GCTAGCTAGCTAGCTA"""
    
    csv_path = os.path.join(temp_dir, "blastn_output.csv")
    with open(csv_path, "w") as f:
        f.write(csv_content)
    
    return csv_path


@pytest.fixture
def sample_csv_blastx(temp_dir):
    """Create a sample BLAST output CSV file for BLASTx testing."""
    csv_content = """query1	gi|111111	98.0	50	1	0	1	50	1	50	-51	185	123456	Archaea	1	1	150	ATGCATGCATGCATGC	ATGCATGCATGCATGC
query1	gi|222222	95.0	50	2	1	1	50	50	100	1e-40	155	654321	Bacteria	1	1	150	ATGCATGCATGCATGC	ATGAATGCATGCATGC"""
    
    csv_path = os.path.join(temp_dir, "blastx_output.csv")
    with open(csv_path, "w") as f:
        f.write(csv_content)
    
    return csv_path


@pytest.fixture
def sample_config_file(temp_dir):
    """Create a sample d-chimer config file for testing."""
    config_content = """blastn_parameters:
  dbpath_nt: /path/to/nt/nt
  nb_threads_bn: 4
  evalue_nt: 0.01

filter_blastn_parameters:
  d: 10
  l: 50

blastx_parameters:
  dbpath_vrl: /path/to/viral/db
  dbpath_nr: /path/to/nr/nr
  nb_threads_bx: 4
  evalue_vir: 0.1
  evalue_nr: 0.01

filter_blastx_parameters:
  d: 10
  l: 17

add_taxo_parameters:
  tax_lineages_file: /path/to/taxonomy/file.dmp

blast_path: /usr/bin
"""
    
    config_path = os.path.join(temp_dir, "dchimer_config.yaml")
    with open(config_path, "w") as f:
        f.write(config_content)
    
    return config_path
