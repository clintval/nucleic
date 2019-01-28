from collections import OrderedDict

__all__ = ['MutRecord']


class MutRecord(OrderedDict):
    """A record from a `.mut` file.

    Note:
        See description of this record at the following URL:

            - https://software.broadinstitute.org/software/igv/MUT
    """

    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    def from_line(line):
        pass

    @staticmethod
    def contig(self):
        return self['Contig']

    @staticmethod
    def start(self):
        return self['Start']

    @staticmethod
    def end(self):
        return self['End']

    @staticmethod
    def sample(self):
        return self['Sample']

    @staticmethod
    def var_type(self):
        return self['VariationType']
    
    @staticmethod
    def ref(self):
        return self['REF']

    @staticmethod
    def alt(self):
        return self['ALT']
