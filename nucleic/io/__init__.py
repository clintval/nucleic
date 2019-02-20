from collections import OrderedDict
from typing import Any, List, Type

__all__ = ['MutRecord']

class MutRecordKeys(object):
    contig = 'contig'
    start = 'start'
    end = 'end'
    sample = 'sample'
    variation_type = 'variation_type'
    ref = 'ref'
    alt = 'alt'
    alt_depth = 'alt_depth'
    depth = 'depth'
    dp_ad_difference = 'dp_ad_difference'
    subtype = 'subtype'
    context = 'context'
    id_field = 'id'
    filter_field = 'filter'
    qual_field = 'qual'


class MutRecord(OrderedDict):
    """A record from a `.mut` file.

    Note:
        See description of this record at the following URL:

            - https://software.broadinstitute.org/software/igv/MUT
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

    @classmethod
    def from_line(self, line: List[str]) -> 'MutRecord':
        return MutRecord()

    @property
    def contig(self: 'MutRecord') -> str:
        """Return the name of the reference sequence."""
        contig: str = self['Contig']
        return contig

    @property
    def start(self: 'MutRecord') -> int:
        """Return the 0-based start coordinate."""
        start: int = self['Start']
        return start

    @property
    def end(self: 'MutRecord') -> int:
        """Return the half-open end coordinate."""
        end: int = self['End']
        return end

    @property
    def sample(self: 'MutRecord') -> str:
        """Return the sample name."""
        sample: str = self['Sample']
        return sample

    @property
    def var_type(self: 'MutRecord') -> str:
        """Return this type of variation."""
        var_type: str = self['VariationType']
        return var_type

    @property
    def ref(self: 'MutRecord') -> str:
        """Return the reference allele."""
        ref: str = self['REF']
        return ref

    @property
    def alt(self: 'MutRecord') -> str:
        """Return the alternate allele."""
        alt: str = self['ALT']
        return alt
