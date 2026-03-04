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


class Subject():
    """
    Represent a single BLAST hit/match with alignment details.

    This class parses and stores information about a single subject sequence match
    from a BLAST output line. Handles both BLASTn (nucleotide) and BLASTx (protein)
    output formats with appropriate field parsing and sense direction detection.

    Attributes:
        sacc: Subject sequence accession number
        pident: Percent identity
        length: Alignment length
        mismatch: Number of mismatches
        gaps: Number of gaps
        qstart: Query sequence start coordinate
        qend: Query sequence end coordinate
        sstart: Subject sequence start coordinate
        send: Subject sequence end coordinate
        evalue: E-value of the match
        bitscore: Bit score of the match
        staxid: Subject taxonomy ID
        sskingdom: Subject kingdom classification
        qlen: Query sequence length
        qseq: Query sequence alignment
        sseq: Subject sequence alignment
        sens: Strand/frame sense direction (True = forward, False = reverse)
        sframe: Subject frame (BLASTx only)
        qframe: Query frame (BLASTx only)
    """

    def __init__(self, bl, program):
        """
        Initialize a Subject from a BLAST output line.

        Parses BLAST output columns according to the program type (BLASTn or BLASTx)
        and determines sense direction for strand/frame orientation.

        Args:
            bl: List of BLAST output fields (tab-separated line split)
            program: BLAST program type ('blastn' or 'blastx')
        """
        if program == 'blastx':
            # For BLASTx, check qframe (column 15) to determine sense direction
            # Negative frame indicates reverse strand
            if int(bl[15]) < 0:
                self.sacc, self.pident, self.length, self.mismatch, \
                    self.gaps, self.qend, self.qstart, self.sstart, \
                    self.send, self.evalue, self.bitscore, self.staxid, \
                    self.sskingdom, self.sframe, self.qframe, self.qlen, \
                    self.qseq, self.sseq = bl[1:19]
                self.sens = False
            else:
                self.sacc, self.pident, self.length, self.mismatch, \
                    self.gaps, self.qstart, self.qend, self.send, \
                    self.sstart, self.evalue, self.bitscore, self.staxid, \
                    self.sskingdom, self.sframe, self.qframe, self.qlen, \
                    self.qseq, self.sseq = bl[1:19]
                self.sens = True

        elif program == 'blastn':
            # For BLASTn, compare sstart/send to determine sense direction
            # If sstart > send, sequence is on reverse strand
            self.sacc, self.pident, self.length, \
                self.mismatch, self.gaps, self.qstart, \
                self.qend, self.sstart, self.send, \
                self.evalue, self.bitscore, self.staxid, \
                self.sskingdom, self.qlen, self.qseq, self.sseq = bl[1:17]
            self.sens = True
            if int(bl[9]) < int(bl[8]):
                self.sens = False
