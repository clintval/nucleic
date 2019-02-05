from collections import OrderedDict
from typing import Any, List, Type

__all__ = ['MutRecord']


class MutRecord(OrderedDict):
    """A record from a `.mut` file.

    Note:
        See description of this record at the following URL:

            - https://software.broadinstitute.org/software/igv/MUT
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        pass

    @classmethod
    def from_line(self, line: List[str]) -> 'MutRecord':
        return MutRecord()

    @staticmethod
    def contig(self: 'MutRecord') -> str:
        """Return the name of the reference sequence."""
        contig: str = self['Contig']
        return contig

    @staticmethod
    def start(self: 'MutRecord') -> int:
        """Return the 0-based start coordinate."""
        start: int = self['Start']
        return start

    @staticmethod
    def end(self: 'MutRecord') -> int:
        """Return the half-open end coordinate."""
        end: int = self['End']
        return end

    @staticmethod
    def sample(self: 'MutRecord') -> str:
        """Return the sample name."""
        sample: str = self['Sample']
        return sample

    @staticmethod
    def var_type(self: 'MutRecord') -> str:
        """Return this type of variation."""
        var_type: str = self['VariationType']
        return var_type

    @staticmethod
    def ref(self: 'MutRecord') -> str:
        """Return the reference allele."""
        ref: str = self['REF']
        return ref

    @staticmethod
    def alt(self: 'MutRecord') -> str:
        """Return the alternate allele."""
        alt: str = self['ALT']
        return alt
