import csv

from pathlib import Path
from typing import Iterable, List, Optional, TextIO, Union

import attr

from pyfaidx import Fasta

from nucleic import Dna
from nucleic.io.fasta import subseq


@attr.s(auto_attribs=True, slots=True)
class Interval(object):
    """An interval which can be placed on a contig within a reference."""

    contig: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    name: Optional[str] = None

    @classmethod
    def fields(cls) -> List[str]:
        """List the fields of this class."""
        fields: List[str] = [attribute.name for attribute in attr.fields(cls)]
        return fields

    @classmethod
    def from_iterable(cls, iterable: Iterable) -> 'Inteval':
        """Build this interval from a list."""
        contig, start, end, name, *_ = iterable
        return cls(str(contig), int(start), int(end), str(name))

    def reference_seq(self, reference: Union[Fasta, Path]) -> Dna:
        """Get thee reference sequence of this interval from a FASTA."""
        assert self.is_closed(), f"Interval must have start and end coordinates: {self}"
        assert self.is_placed(), f"Interval have `contig` defined: {self}"
        return subseq(reference, self.contig, self.start, self.end)
    
    def is_closed(self):
        """Test if this is a closed interval."""
        return all([num is not None for num in (self.start, self.end)])

    def is_placed(self):
        """Test if this interval is located on a contig."""
        return self.contig is not None


class IntervalList(list):
    """A container of intervals."""

    def __init__(self, *args):
        super().__init__(args)

    @classmethod
    def read_bed(cls, handle: TextIO):
        """Read a BED filehandle into an interval list."""
        reader = csv.reader(handle, delimiter='\t')
        return cls(*(Interval.from_iterable(line) for line in reader))

    def __repr__(self):
        items = ', '.join([str(item) for item in self])
        return f'{self.__class__.__qualname__}({items})'