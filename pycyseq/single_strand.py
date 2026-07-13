from dataclasses import dataclass

import pysam

from pycyseq.alignment import Alignment


@dataclass
class SingleStrandAlignment:
    """Alignment of a single-strand DNA molecule.

    Args:
        aln: Alignment object
    """

    def __init__(self, aln: Alignment):
        self.aln = aln

    def __str__(self) -> str:
        return self.aln.__str__()

    @property
    def pysam(self) -> pysam.AlignedSegment:
        """The pysam.AlignmentSegment object for access to its methods. See
        https://pysam.readthedocs.io/en/stable/index.html for documentation.
        """
        return self.aln.pysam

    @property
    def name(self) -> str | None:
        """Name of the read"""
        return self.aln._pysam.query_name

    @property
    def has_rotations(self) -> bool:
        """Whether the DNA sequence has more than on possible rotation based on
        alignment.
        """
        return True

    @property
    def num_rotations(self) -> int:
        """The number of possible rotations."""
        return self.aln._pysam.get_tag("YR")  # ty: ignore

    @property
    def mapq(self) -> int:
        """Mapping quality of the alignment."""
        return self.aln._pysam.mapping_quality

    @property
    def len(self) -> int:
        """Number of aligned bases in the DNA sequence"""
        return self.pysam.reference_end - self.pysam.reference_start  # ty: ignore

    @property
    def is_duplicate(self) -> int:
        """Return true if the alignment has been marked as duplicate"""
        return self.pysam.is_duplicate
