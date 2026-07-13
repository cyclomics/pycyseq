from dataclasses import dataclass

import pysam


@dataclass
class Alignment:
    """Wrapper class around a pysam alignment.

    Args:
        pysam: `pysam.AlignedSegment` object with the alignment information.
    """

    def __init__(self, pysam: pysam.AlignedSegment):
        self._pysam = pysam

    @property
    def pysam(self) -> pysam.AlignedSegment:
        """The pysam.AlignmentSegment object for access to its methods. See
        https://pysam.readthedocs.io/en/stable/index.html for documentation.
        """
        return self._pysam

    def __str__(self) -> str:
        """Print the alignment coordinates and direction."""

        if self._pysam.is_mapped:
            return f"{self._pysam.reference_name}-{self._pysam.reference_start}:{self._pysam.reference_end}-{self.direction}"
        else:
            return "Not mapped"

    @property
    def direction(self) -> str:
        """Direction of the alignment of the read, can be 'fwd', 'rev' or 'none' if not aligned."""
        if self._pysam.is_mapped:
            if self._pysam.is_forward:
                return "fwd"
            else:
                return "rev"
        else:
            return "none"

    @property
    def is_mapped(self) -> bool:
        """Return true if the alignment is mapped."""
        return self._pysam.is_mapped
