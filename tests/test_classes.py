"""Unit tests for d-chimer core classes (Subject, Contig, CsvIO)."""

import pytest
from io import StringIO
from dchimer.Subject import Subject
from dchimer.Contig import Contig
from dchimer.CsvIO import CsvIO


class TestSubject:
    """Tests for Subject class."""

    def test_blastn_subject_forward_strand(self):
        """Test Subject initialization from BLASTn output (forward strand)."""
        # BLASTn output line: qseqid, sacc, pident, length, mismatch, gaps,
        # qstart, qend, sstart, send, evalue, bitscore, staxid, sskingdom, qlen, qseq, sseq
        bl = [
            "query1",  # qseqid (index 0, skipped)
            "gi|123456",  # sacc (1)
            "95.5",  # pident (2)
            "100",  # length (3)
            "4",  # mismatch (4)
            "1",  # gaps (5)
            "1",  # qstart (6)
            "100",  # qend (7)
            "1",  # sstart (8)
            "100",  # send (9)
            "1.5e-50",  # evalue (10)
            "185",  # bitscore (11)
            "562",  # staxid (12)
            "Bacteria",  # sskingdom (13)
            "100",  # qlen (14)
            "ATGC" * 25,  # qseq (15)
            "ATGC" * 25,  # sseq (16)
        ]
        
        subject = Subject(bl, 'blastn')
        assert subject.sacc == "gi|123456"
        assert subject.pident == "95.5"
        assert subject.length == "100"
        assert subject.sens == True  # Forward strand (sstart < send)
        assert subject.staxid == "562"

    def test_blastn_subject_reverse_strand(self):
        """Test Subject initialization from BLASTn output (reverse strand)."""
        bl = [
            "query1",
            "gi|789012",
            "92.0",
            "100",
            "8",
            "0",
            "1",
            "100",
            "100",  # sstart > send = reverse strand
            "1",    # send
            "3.2e-40",
            "155",
            "562",
            "Bacteria",
            "100",
            "ATGC" * 25,
            "ATGC" * 25,
        ]
        
        subject = Subject(bl, 'blastn')
        assert subject.sens == False  # Reverse strand

    def test_blastx_subject_forward_frame(self):
        """Test Subject initialization from BLASTx output (positive frame)."""
        # BLASTx has additional fields: sframe, qframe
        bl = [
            "query1",  # index 0, skipped
            "gi|111111",  # sacc (1)
            "98.0",  # pident (2)
            "50",  # length (3)
            "1",  # mismatch (4)
            "0",  # gaps (5)
            "1",  # qstart (6)
            "50",  # qend (7)
            "1",  # send (8) 
            "50",  # sstart (9)
            "-51",  # evalue (10)
            "185",  # bitscore (11)
            "123456",  # staxid (12)
            "Archaea",  # sskingdom (13)
            "1",  # sframe (14)
            "1",  # qframe (15) - positive = forward
            "150",  # qlen (16)
            "ATGC" * 6,  # qseq (17)
            "ATGC" * 6,  # sseq (18)
        ]
        
        subject = Subject(bl, 'blastx')
        assert subject.sens == True  # Positive frame
        assert subject.qframe == "1"
        assert subject.sframe == "1"

    def test_blastx_subject_reverse_frame(self):
        """Test Subject initialization from BLASTx output (negative frame)."""
        bl = [
            "query1",
            "gi|222222",
            "95.0",
            "50",
            "2",
            "1",
            "1",
            "50",
            "50",
            "100",
            "1e-40",
            "155",
            "654321",
            "Bacteria",
            "-1",  # sframe
            "-1",  # qframe - negative = reverse
            "150",
            "ATGC" * 6,
            "ATGC" * 6,
        ]
        
        subject = Subject(bl, 'blastx')
        assert subject.sens == False  # Negative frame


class TestContig:
    """Tests for Contig class."""

    def test_contig_initialization(self):
        """Test Contig object initialization."""
        contig = Contig("query1", minDec=10, minBpLength=50)
        assert contig.id == "query1"
        assert contig.minDec == 10
        assert contig.minBpLength == 50
        assert len(contig.subject) == 0
        assert len(contig.selectedSubjects) == 0

    def test_contig_attributes(self):
        """Test that Contig attributes are properly typed."""
        contig = Contig("seq123", minDec="20", minBpLength="100")
        assert isinstance(contig.minDec, int)
        assert isinstance(contig.minBpLength, int)
        assert contig.minDec == 20
        assert contig.minBpLength == 100


class TestCsvIO:
    """Tests for CsvIO class."""

    def test_csvio_initialization(self):
        """Test CsvIO object initialization."""
        handle = StringIO("")
        csv_io = CsvIO(handle, minDec=10, minBpLength=50, program='blastn')
        assert csv_io.minDec == 10
        assert csv_io.minBpLength == 50
        assert csv_io.program == 'blastn'

    def test_csvio_next_no_match(self):
        """Test CsvIO.next() when contig is not in output."""
        csv_content = "query1\tgi|123\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t185\t562\tBacteria\t100\tATGC\tATGC\n"
        handle = StringIO(csv_content)
        csv_io = CsvIO(handle, minDec=10, minBpLength=50, program='blastn')
        
        # Looking for query2 which doesn't exist
        result = csv_io.next("query2")
        assert result is None

    def test_csvio_next_empty_file(self):
        """Test CsvIO.next() on empty file."""
        handle = StringIO("")
        csv_io = CsvIO(handle, minDec=10, minBpLength=50, program='blastn')
        
        result = csv_io.next("query1")
        assert result is None

    def test_csvio_with_invalid_data_types(self):
        """Test CsvIO with string parameters that should be integers."""
        handle = StringIO("")
        # String values should be converted to int
        csv_io = CsvIO(handle, minDec="15", minBpLength="75", program='blastn')
        assert csv_io.minDec == 15
        assert csv_io.minBpLength == 75
