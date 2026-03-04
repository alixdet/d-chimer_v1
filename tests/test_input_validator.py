"""Unit tests for dchimer_input_validator module."""

import pytest
import os
import tempfile
from dchimer.dchimer_input_validator import (
    validate_fasta_file,
    validate_fasta_sequences,
)


class TestValidateFastaFile:
    """Tests for validate_fasta_file function."""

    def test_valid_fasta_file(self, sample_fasta_file):
        """Test validation of a valid FASTA file."""
        num_seqs, file_size = validate_fasta_file(sample_fasta_file)
        assert num_seqs == 3
        assert file_size > 0

    def test_nonexistent_file(self):
        """Test error handling for non-existent file."""
        with pytest.raises(FileNotFoundError):
            validate_fasta_file("/nonexistent/path/file.fasta")

    def test_empty_file(self, temp_dir):
        """Test error handling for empty file."""
        empty_file = os.path.join(temp_dir, "empty.fasta")
        open(empty_file, "w").close()
        
        with pytest.raises(ValueError, match="empty"):
            validate_fasta_file(empty_file)

    def test_invalid_fasta_format(self, temp_dir):
        """Test error handling for invalid FASTA format."""
        invalid_fasta = os.path.join(temp_dir, "invalid.fasta")
        with open(invalid_fasta, "w") as f:
            f.write("not a fasta file\njust some text\n")
        
        # Empty content results in ValueError
        with pytest.raises(ValueError, match="No sequences"):
            validate_fasta_file(invalid_fasta)

    def test_sequence_with_empty_sequence(self, temp_dir):
        """Test error handling for sequences with no actual sequence data."""
        fasta_file = os.path.join(temp_dir, "empty_seq.fasta")
        with open(fasta_file, "w") as f:
            f.write(">seq1\n\n>seq2\nATGC\n")
        
        with pytest.raises(ValueError, match="empty"):
            validate_fasta_file(fasta_file)


class TestValidateFastaSequences:
    """Tests for validate_fasta_sequences function."""

    def test_blastn_valid_sequences(self, sample_fasta_file):
        """Test validation of valid nucleotide sequences for BLASTn."""
        stats = validate_fasta_sequences(sample_fasta_file, program='blastn')
        assert stats['num_sequences'] == 3
        assert stats['total_length'] == 34  # 16 + 16 + 2
        assert stats['min_length'] == 2
        assert stats['max_length'] == 16
        assert not stats['has_invalid_chars']

    def test_blastx_amino_acid_sequences(self, temp_dir):
        """Test validation of amino acid sequences for BLASTx."""
        fasta_file = os.path.join(temp_dir, "protein.fasta")
        with open(fasta_file, "w") as f:
            f.write(">prot1\nMVKLLAS\n>prot2\nFEQWER\n")
        
        stats = validate_fasta_sequences(fasta_file, program='blastx')
        assert stats['num_sequences'] == 2
        assert stats['total_length'] == 13
        assert not stats['has_invalid_chars']

    def test_invalid_nucleotide_chars(self, temp_dir):
        """Test error handling for invalid nucleotide characters."""
        fasta_file = os.path.join(temp_dir, "invalid_nuc.fasta")
        with open(fasta_file, "w") as f:
            f.write(">seq1\nATGCXYZ\n")  # X, Y, Z are invalid for nucleotides
        
        with pytest.raises(ValueError, match="Invalid character"):
            validate_fasta_sequences(fasta_file, program='blastn')

    def test_invalid_amino_acid_chars(self, temp_dir):
        """Test error handling for invalid amino acid characters."""
        fasta_file = os.path.join(temp_dir, "invalid_aa.fasta")
        with open(fasta_file, "w") as f:
            f.write(">prot1\nMVK2LAS\n")  # '2' is not a valid amino acid or nucleotide
        
        with pytest.raises(ValueError, match="Invalid character"):
            validate_fasta_sequences(fasta_file, program='blastx')

    def test_valid_iupac_nucleotide_codes(self, temp_dir):
        """Test that IUPAC ambiguity codes are accepted for nucleotides."""
        fasta_file = os.path.join(temp_dir, "iupac.fasta")
        with open(fasta_file, "w") as f:
            # N = any base, W = A or T, S = G or C, etc.
            f.write(">seq1\nATGCNNWWSS\n")
        
        stats = validate_fasta_sequences(fasta_file, program='blastn')
        assert stats['num_sequences'] == 1
        assert not stats['has_invalid_chars']

    def test_stop_codon_allowed_in_proteins(self, temp_dir):
        """Test that stop codon (*) is accepted in protein sequences."""
        fasta_file = os.path.join(temp_dir, "with_stop.fasta")
        with open(fasta_file, "w") as f:
            f.write(">prot1\nMVKLLAS*\n")  # * = stop codon
        
        stats = validate_fasta_sequences(fasta_file, program='blastx')
        assert stats['num_sequences'] == 1
        assert not stats['has_invalid_chars']
