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

from . import Contig
from . import Subject


class CsvIO():
    """
    Parse BLAST CSV output and return contigs with their subject matches.

    This class reads BLAST CSV output line by line and groups output by contig ID,
    creating Contig objects and associating Subject objects with each match.

    Attributes:
        handle: File handle to BLAST CSV output
        minDec: Minimum number of hits required (filter parameter)
        minBpLength: Minimum alignment length required (filter parameter)
        program: BLAST program type ('blastn' or 'blastx')
    """

    def __init__(self, handle, minDec, minBpLength, program):
        """
        Initialize a BLAST CSV output parser.

        Args:
            handle: File handle to BLAST CSV output
            minDec: Minimum number of hits to consider (int)
            minBpLength: Minimum alignment length in base pairs (int)
            program: BLAST program type ('blastn' or 'blastx')
        """
        self.handle = handle
        self.minDec = int(minDec)
        self.minBpLength = int(minBpLength)
        self.program = program

    def next(self, s_id):
        """
        Get next contig from BLAST output and its matches.

        Reads lines from BLAST CSV output until all matches for a given contig ID
        are collected, then returns a Contig object with all its Subject matches.

        Args:
            s_id: Contig ID to search for

        Returns:
            Contig object with matches if found, None if contig not in output
        """
        pos = self.handle.tell()
        bl = self.handle.readline().split()

        # No BLAST output line (end of file)
        if not bl:
            return None

        c = Contig.Contig(bl[0], self.minDec, self.minBpLength)  # init Contig
        # if Current contig (from fasta), does not appear in (output) csv file
        if c.id != s_id:
            self.handle.seek(pos)  # Return to initial position
            return None

        c.subject.append(Subject.Subject(bl, self.program))

        # Collect subsequent lines with the same query contig id
        while c.id == s_id:
            pos = self.handle.tell()
            bl = self.handle.readline().split()
            if not bl:
                return c  # end of file
            c.id = bl[0]
            c.subject.append(Subject.Subject(bl, self.program))

        # End of this contig's matches
        c.id = s_id
        c.subject = c.subject[:-1]
        self.handle.seek(pos)

        return c
