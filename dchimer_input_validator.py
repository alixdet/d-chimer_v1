#!/usr/bin/python3
# -*-coding:utf-8 -*

"""
Input validation for d-chimer.

Validates FASTA input files before processing.
"""

import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaParseError


def validate_fasta_file(fasta_path):
    """
    Validate that a FASTA file exists, is readable, and contains valid sequences.

    Args:
        fasta_path: path to FASTA file

    Returns:
        tuple: (num_sequences, file_size_kb)

    Raises:
        FileNotFoundError: if file doesn't exist or is empty
        ValueError: if file is not valid FASTA format
    """
    # Check if file exists
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_path}")

    # Check if file is readable
    if not os.access(fasta_path, os.R_OK):
        raise PermissionError(f"Cannot read FASTA file: {fasta_path}")

    # Check if file is empty
    file_size = os.path.getsize(fasta_path)
    if file_size == 0:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    # Try to parse FASTA file to validate format
    try:
        num_sequences = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            if not record.seq:
                raise ValueError(
                    f"Sequence '{record.id}' in {fasta_path} is empty"
                )
            num_sequences += 1

        if num_sequences == 0:
            raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

        return num_sequences, file_size / 1024  # Return size in KB

    except FastaParseError as e:
        raise ValueError(
            f"Invalid FASTA format in {fasta_path}: {e}"
        )
    except Exception as e:
        raise ValueError(
            f"Error parsing FASTA file {fasta_path}: {e}"
        )


def validate_fasta_sequences(fasta_path, program='blastn'):
    """
    Additional validation specific to BLAST program type.

    Args:
        fasta_path: path to FASTA file
        program: BLAST program type ('blastn' or 'blastx')

    Returns:
        dict: statistics about sequences

    Raises:
        ValueError: if sequences are invalid for the program type
    """
    stats = {
        'num_sequences': 0,
        'total_length': 0,
        'min_length': float('inf'),
        'max_length': 0,
        'has_invalid_chars': False,
    }

    # Valid nucleotides (BLASTN)
    valid_nucleotides = set('ATGCNatgcnWwSsMmKkRrYyBbDdHhVv')
    # Valid amino acids (BLASTX)
    valid_amino_acids = set('ACDEFGHIKLMNPQSTVWY*acdefghiklmnpqstvwy')

    valid_chars = valid_nucleotides if program == 'blastn' else valid_amino_acids

    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_len = len(record.seq)
            stats['num_sequences'] += 1
            stats['total_length'] += seq_len
            stats['min_length'] = min(stats['min_length'], seq_len)
            stats['max_length'] = max(stats['max_length'], seq_len)

            # Check for invalid characters
            for char in str(record.seq):
                if char not in valid_chars:
                    stats['has_invalid_chars'] = True
                    raise ValueError(
                        f"Invalid character '{char}' in sequence '{record.id}' "
                        f"for {program} (file: {fasta_path})"
                    )

        if stats['min_length'] == float('inf'):
            stats['min_length'] = 0

        return stats

    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Error validating sequences in {fasta_path}: {e}")
