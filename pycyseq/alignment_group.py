from typing import List
from dataclasses import dataclass

from pycyseq import SingleStrandAlignment, DoubleStrandAlignment


@dataclass
class AlignmentGroup:
    """Group of alignments that belong to the same CySeq read.

    Args:
        alignments: List of SingleStrandAlignment or DoubleStrandAlignment objects.

    Raises:
        ValueError: If not all objects in alignments are of the same type
            SingleStrandAlignment or DoubleStrandAlignment.
    """

    def __init__(
        self, alignments: List[SingleStrandAlignment] | List[DoubleStrandAlignment]
    ):
        if len(alignments) == 0:
            raise ValueError

        ssdna = 0
        dsdna = 0
        for aln in alignments:
            if isinstance(aln, SingleStrandAlignment):
                ssdna += 1
                continue
            if isinstance(aln, DoubleStrandAlignment):
                dsdna += 1
                continue
            raise ValueError

        if ssdna > 0 and dsdna > 0:
            raise ValueError

        self.alignments = alignments

    def __str__(self) -> str:
        if self.is_ssdna:
            return f"Single-strand group with {len(self)} alignment(s)"

        if self.is_dsdna:
            return f"Double-strand group with {len(self)} alignment(s)"

        return "Unknown"

    def __len__(self) -> int:
        """Number of alignments."""
        return len(self.alignments)

    def __iter__(self):
        return iter(self.alignments)

    def __getitem__(self, index):
        return self.alignments[index]

    @property
    def name(self) -> str | None:
        """Name of the CySeq read."""
        return self.alignments[0].name

    @property
    def is_ssdna(self) -> bool:
        """Return true if the alignments in this object are from ssDNA molecules."""
        if isinstance(self.alignments[0], SingleStrandAlignment):
            return True
        else:
            return False

    @property
    def is_dsdna(self) -> bool:
        """Return true if the alignments in this object are from dsDNA molecules."""
        return not self.is_ssdna

    @property
    def has_rotations(self) -> bool:
        """Return true if there is more than one possible alignment."""
        return len(self.alignments) > 1
