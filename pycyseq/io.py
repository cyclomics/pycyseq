from typing import List, Iterator
from pathlib import Path

import pysam

from pycyseq import (
    SingleStrandAlignment,
    DoubleStrandAlignment,
    AlignmentGroup,
    Alignment,
)


def wrap_aligments(alns: List[Alignment]) -> AlignmentGroup:
    read_types = set()
    for aln in alns:
        read_types.add(aln.pysam.get_tag("XT"))

    if len(read_types) > 2:
        raise ValueError(f"Mix of read types: {read_types}")

    ssdna = 0
    dsdna = 0
    for rt in read_types:
        if "1D" in rt or "2D" in rt:
            ssdna += 1
        elif "3D" in rt or "4D" in rt:
            dsdna += 1
        else:
            continue

    if dsdna >= 1 and ssdna >= 1:
        raise ValueError

    if dsdna == 0 and ssdna == 0:
        raise ValueError

    if ssdna >= 1 and dsdna == 0:
        new_alns: List[SingleStrandAlignment] = []
        for aln in alns:
            new_alns.append(SingleStrandAlignment(aln))

    if dsdna >= 1 and ssdna == 0:
        num_rotations = len(alns) // 2
        fwd = [None] * num_rotations
        rev = [None] * num_rotations

        for aln in alns:
            rotation_number: int = aln.pysam.get_tag("YO")  # ty: ignore
            if aln.pysam.is_forward:
                fwd[rotation_number] = aln
            else:
                rev[rotation_number] = aln

        new_alns: List[DoubleStrandAlignment] = []
        for f, r in zip(fwd, rev):
            new_alns.append(DoubleStrandAlignment(fwd=f, rev=r))

    return AlignmentGroup(new_alns)


def read_cyseq_alignment_file(
    file: Path | str, only_successful: bool = True
) -> Iterator[AlignmentGroup]:
    """Read the data from a CySeq sam or bam file.

    Parameters:
        file: sam or bam file to read
        only_successful: skip reads that failed CySeq quality control. For example,
            consensus was not generated, or did not align to the reference genome.

    Raises:
        ValueError: if the given file does not end in `.sam` or `.bam`.
        NotImplementedError: if only_successful is false.

    Yields:
        An `AlignmentGroup` object with all alignments related to the same
            CySeq read.
    """

    if not only_successful:
        raise NotImplementedError(
            "Retrieval of non-aligned are failed reads not implemented"
        )

    if isinstance(file, str):
        file = Path(file)

    if file.suffix == ".sam":
        mode = "r"
    elif file.suffix == ".bam":
        mode = "rb"
    else:
        raise ValueError(
            f"Unexpected file extension. Expected '.sam' or '.bam'. Received: {file.suffix}"
        )

    with pysam.AlignmentFile(str(file), mode) as f:
        alns = []
        current_read_id = None

        for aln in f:
            alignment = Alignment(pysam=aln)
            if current_read_id is None:
                current_read_id = alignment.pysam.query_name

            # New read encountered → flush previous
            if alignment.pysam.query_name != current_read_id:
                if len(alns) > 0:
                    yield wrap_aligments(alns)

                alns = []
                current_read_id = alignment.pysam.query_name

            if only_successful:
                if alignment.pysam.is_qcfail:
                    continue
                else:
                    alns.append(alignment)
            else:
                alns.append(alignment)

        if len(alns) > 0:
            yield wrap_aligments(alns)
