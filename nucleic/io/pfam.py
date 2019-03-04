import json
import urllib.request as request
from typing import Dict, List, Mapping, Optional, Tuple

import attr

__all__ = ['PfamGraphicFeature', 'PfamGraphicResponse']


@attr.s(auto_attribs=True, slots=True)
class PfamGraphicFeature(object):
    """An object representing a Pfam graphic.

    Args:
        feature: A dictionary of graphic feature options.

    """
    color: str = 'grey'
    display: Optional[str] = None
    end: Optional[int] = None
    endstyle: Optional[str] = None
    link: Optional[str] = None
    start: Optional[int] = None
    startstyle: Optional[str] = None
    text: Optional[str] = ''
    feature_type: Optional[str] = None
    metadata: Dict = dict()


class PfamGraphicResponse(object):
    """An object representing a Pfam graphic server response."""

    def __init__(self, content: Mapping) -> None:
        self.regions: List[PfamGraphicFeature] = []

        self.length: Optional[int] = content.get('length', None)
        self.markups: Tuple[str] = content.get('markups', tuple())
        self.metadata: Dict = content.get('metadata', dict())
        self.motifs: Tuple[str] = content.get('motifs', tuple())

        for feature in content.get('regions', tuple()):
            self.regions.append(PfamGraphicFeature(**feature))

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}(\n'
            f'  accession:    "{self.metadata.get("accession", "N/A")}"\n'
            f'  identifier:   "{self.metadata.get("identifier", "N/A")}"\n'
            f'  organism:     "{self.metadata.get("organism", "N/A")}")\n'
            f'  description:  "{self.metadata.get("description", "N/A")}"\n'
            f'  number of motifs:  {len(self.motifs)}\n'
            f'  number of regions: {len(self.regions)})'
        )


def fetch_pfam_graphics(pfam_id: str) -> List[PfamGraphicResponse]:
    url = f'http://pfam.xfam.org/protein/{pfam_id}/graphic'
    with request.urlopen(url) as handle:
        responses = list(map(PfamGraphicResponse, json.load(handle)))
    return responses
