import io
from collections import OrderedDict
from itertools import groupby
from operator import attrgetter
from pathlib import Path
from typing import Any, Dict, List, Union

from nucleic import Dna, Snv, Variant
from nucleic.metrics import Metric

__all__ = ['MutSnvMetric', 'MutRecord', 'MutReader']


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
    normalized_subtype = 'normalized_subtype'
    normalized_context = 'normalized_context'
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

    @property
    def contig(self: 'MutRecord') -> str:
        """Return the name of the reference sequence."""
        contig: str = self[MutRecordKeys.contig]
        return contig

    @property
    def start(self: 'MutRecord') -> int:
        """Return the 0-based start coordinate."""
        start: str = self[MutRecordKeys.start]
        return int(start)

    @property
    def end(self: 'MutRecord') -> int:
        """Return the half-open end coordinate."""
        end: str = self[MutRecordKeys.end]
        return int(end)

    @property
    def sample(self: 'MutRecord') -> str:
        """Return the sample name."""
        sample: str = self[MutRecordKeys.sample]
        return sample

    @property
    def variation_type(self: 'MutRecord') -> str:
        """Return this type of variation."""
        variation_type: str = self[MutRecordKeys.variation_type]
        return variation_type

    @property
    def ref(self: 'MutRecord') -> Dna:
        """Return the reference allele."""
        ref: str = self[MutRecordKeys.ref]
        return Dna(ref)

    @property
    def alt(self: 'MutRecord') -> Dna:
        """Return the alternate allele."""
        alt: str = self[MutRecordKeys.alt]
        return Dna(alt)

    @property
    def alt_depth(self: 'MutRecord') -> int:
        """Return the alternate allele depth."""
        alt_depth: str = self[MutRecordKeys.alt_depth]
        return int(alt_depth)

    @property
    def depth(self: 'MutRecord') -> int:
        """Return the total allele depth."""
        depth: str = self[MutRecordKeys.depth]
        return int(depth)

    @property
    def dp_ad_difference(self: 'MutRecord') -> int:
        """Return the difference of all alternate allele depths from the total depth."""
        dp_ad_difference: str = self[MutRecordKeys.dp_ad_difference]
        return int(dp_ad_difference)

    @property
    def normalized_subtype(self: 'MutRecord') -> str:
        """Return the normalized subtitution type."""
        normalized_subtype: str = self[MutRecordKeys.normalized_subtype]
        return normalized_subtype

    @property
    def normalized_context(self: 'MutRecord') -> Dna:
        """Return the normalized local context."""
        normalized_context: str = self[MutRecordKeys.normalized_context]
        return Dna(normalized_context)

    def to_variant(self: 'MutRecord') -> Union[Variant, Snv]:
        """Return this record as a :class:`Variant`."""
        variant = Variant(ref=self.ref, alt=self.alt, data=self)
        if isinstance(variant, Snv):
            variant.context = self.normalized_context
        return variant


class MutReader(object):
    """Read a MUT file."""

    attributes = [
        'contig',
        'start',
        'end',
        'sample',
        'variation_type',
        'ref',
        'alt',
        'alt_depth',
        'depth',
        'dp_ad_difference',
        'normalized_subtype',
        'normalized_context',
    ]

    def __init__(self, handle: io.IOBase) -> None:
        self.handle = handle
        self.header = next(handle).strip().split()
        if not self.header[: len(self.attributes)] == self.attributes:
            raise ValueError('Required column names do not exist: {", ".join(self.attributes)}')

    def __iter__(self) -> 'MutReader':
        return self

    def __next__(self) -> MutRecord:
        """Iterate through rows as dictionaries, then unpack into a :class:`MutRecord`."""
        fields = next(self.handle).strip().split()
        mapping = dict(zip(self.header, fields))
        return MutRecord(**mapping)  # type: ignore

    @staticmethod
    def read_from_path(path: Path) -> Dict[str, List[MutRecord]]:
        """Create a map of samples to :class:`MutRecord` from a file path."""
        with open(path, 'r') as handle:
            sample_map = {
                name: list(group)
                for name, group in groupby(MutReader(handle), attrgetter('sample'))  # type: ignore
            }
        return sample_map

    def __repr__(self) -> str:
        return f'{self.__class__.__qualname__}({self.handle})'


class MutSnvMetric(Metric):
    """A metric class holding information about single nucleotide variants."""

    def __init__(self, num_records: int = 0, total_alt_depth: int = 0) -> None:
        self.num_records = num_records
        self.total_alt_depth = total_alt_depth

    def update(self, record: MutRecord) -> None:
        """Update this metric with a new :class:`MutRecord`."""
        if not record.to_variant().is_snv():
            return None
        self.num_records += 1
        self.total_alt_depth += record.alt_depth
        return None
