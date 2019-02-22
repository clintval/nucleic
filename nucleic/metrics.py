import inspect
from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, List


class Metric(ABC):
    """Abstract base class for any metric containing class."""

    @abstractmethod
    def __init__(self) -> None:
        return None

    @abstractmethod
    def update(self, item: Any) -> None:
        """Update this metric from an item."""
        return None

    def update_from(self, iterable: Iterable[Any]) -> None:
        """Update this metric from a collection of items."""
        for item in iterable:
            self.update(item)
        return None

    def fields(self) -> List[str]:
        """List the fields of this class."""
        signature = inspect.signature(self.__init__)  # type: ignore
        fields = sorted(signature.parameters.keys())
        return fields

    def to_dict(self) -> Dict[str, Any]:
        """Transform this class to a dictionaary."""
        mapping = dict(zip(self.fields(), self.to_line()))
        return mapping

    def to_line(self) -> List[Any]:
        """Write the values of this metric to a line."""
        line = [getattr(self, field) for field in self.fields()]
        return line

    def __repr__(self) -> str:
        groups = [f'{field}={repr(key)}' for field, key in self.to_dict().items()]
        return f'{self.__class__.__qualname__}({", ".join(groups)})'

    def __str__(self) -> str:
        return repr(self)


class MetricWriter(object):
    """Not implemented yet."""

    pass
