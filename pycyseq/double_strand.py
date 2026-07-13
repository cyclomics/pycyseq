from dataclasses import dataclass
from typing import Tuple
from functools import cached_property

from pycyseq.alignment import Alignment


@dataclass
class DoubleStrandAlignment:
    """Alignment of a double-strand DNA molecule.

    Args:
        fwd: Alignment object referring to the forward aligning strand
        rev: Alignment object referring to the reverse aligning strand

    Raises:
        ValueError: if the input fails validation
    """

    def __init__(self, fwd: Alignment, rev: Alignment):
        if not isinstance(fwd, Alignment):
            raise ValueError(f"fwd is not an Alignment object, given: {type(fwd)}")

        if not fwd.pysam.is_forward:
            raise ValueError("Alignment given to fwd is reverse aligned")

        if not isinstance(rev, Alignment):
            raise ValueError(f"rev is not an Alignment object, given: {type(rev)}")

        if rev.pysam.is_forward:
            raise ValueError("Alignment given to rev is forward aligned")

        self._fwd = fwd
        self._rev = rev

    def __str__(self) -> str:
        return self._fwd.aln.__str__()  # ty: ignore

    @property
    def fwd(self) -> Alignment:
        """Access the Alignment object of the forward mapping strand."""
        return self._fwd

    @property
    def rev(self) -> Alignment:
        """Access the Alignment object of the reverse mapping strand."""
        return self._rev

    @cached_property
    def type(self) -> str:
        """
        Type of dsDNA molecule

        We classify the dsDNA molecule into nine different types based on the location
        of the overhangs.

        Types of dsDNA:
            - **blunt_double**: there are no overhangs on either end.
            - **blunt_fiveprime_fwd_overhang**: the reference 5' end is blunt,
                with a forward strand overhang on the reference 3'
            - **blunt_fiveprime_rev_overhang**: the reference 5' end is blunt,
                with a reverse strand overhang on the reference 3'
            - **blunt_threeprime_fwd_overhang**: the reference 3' end is blunt,
                with a forward strand overhang on the reference 5'
            - **blunt_threeprime_rev_overhang**: the reference 3' end is blunt,
                with a reverse strand overhang on the reference 5'
            - **jagged_shifted_fiveprime**: there are overhangs on both strand 5' ends
            - **jagged_shifted_threeprime**: there are overhangs on both strand 3' ends
            - **jagged_encompassing_fwd_overhang**: both overhangs are on the forward strand
            - **jagged_encompassing_rev_overhang**: both overhangs are on the reverse strand

        Note that the naming, we use reference 5'/3' and strand 5'/3' to refer to the
            ends relative to the reference or to the strand.

        ??? info "Click to expand for visual representations of each type"

            ```
            Ref 5' -> 3'
            Fwd 5' -> 3'
            Rev 3' <- 5'

            blunt_double
                NNNNNNNNNNNN
                NNNNNNNNNNNN

            blunt_fiveprime_fwd_overhang
                NNNNNNNNNNNN
                NNNNNNNNN___

            blunt_fiveprime_rev_overhang
                NNNNNNNNN___
                NNNNNNNNNNNN

            blunt_threeprime_fwd_overhang
                NNNNNNNNNNNN
                ___NNNNNNNNN

            blunt_threeprime_rev_overhang
                ___NNNNNNNNN
                NNNNNNNNNNNN

            jagged_shifted_fiveprime
                NNNNNNNNN___
                ___NNNNNNNNN

            jagged_shifted_threeprime
                ___NNNNNNNNN
                NNNNNNNNN___

            jagged_encompassing_fwd_overhang
                NNNNNNNNNNNN
                ___NNNNNN___

            jagged_encompassing_rev_overhang
                ___NNNNNN___
                NNNNNNNNNNNN

            ```
        """

        fwd_st = self.fwd.pysam.reference_start
        fwd_en = self.fwd.pysam.reference_end
        rev_st = self.rev.pysam.reference_start
        rev_en = self.rev.pysam.reference_end

        if fwd_st == rev_st:
            if fwd_en == rev_en:
                return "blunt_double"
            elif fwd_en > rev_en:  # ty: ignore
                return "blunt_fiveprime_fwd_overhang"
            else:
                return "blunt_fiveprime_rev_overhang"

        if fwd_en == rev_en:
            if fwd_st < rev_st:
                return "blunt_threeprime_fwd_overhang"
            else:
                return "blunt_threeprime_rev_overhang"

        if fwd_st < rev_st:
            if fwd_en < rev_en:  # ty: ignore
                return "jagged_shifted_fiveprime"
            else:
                return "jagged_encompassing_fwd_overhang"

        if fwd_st > rev_st:
            if fwd_en > rev_en:  # ty: ignore
                return "jagged_shifted_threeprime"
            else:
                return "jagged_encompassing_rev_overhang"

        raise ValueError(f"{fwd_st}-{fwd_en}, {rev_st}-{rev_en}")

    @property
    def length_fwd(self) -> int:
        """Length of the forward mapping strand."""
        fwd_st = self.fwd.pysam.reference_start
        fwd_en = self.fwd.pysam.reference_end
        return fwd_en - fwd_st  # ty: ignore

    @property
    def length_rev(self) -> int:
        """Length of the reverse mapping strand."""
        rev_st = self.rev.pysam.reference_start
        rev_en = self.rev.pysam.reference_end
        return rev_en - rev_st  # ty: ignore

    @property
    def length_sum(self) -> int:
        """Length of the forward and reverse mapping strands summed."""
        return self.length_fwd + self.length_rev

    @property
    def length_long(self) -> int:
        """Length of the longest of the two mapping strands."""
        return max(self.length_fwd, self.length_rev)

    @property
    def length_short(self) -> int:
        """Length of the shortest of the two mapping strands."""
        return min(self.length_fwd, self.length_rev)

    @property
    def longest_strand(self) -> str:
        """Return the direction of the longest strand.

        Can be 'fwd', 'rev', or 'both' if both strands are equally long.
        """
        if self.length_fwd > self.length_rev:
            return "fwd"
        elif self.length_fwd < self.length_rev:
            return "rev"
        else:
            return "both"

    @property
    def shortest_strand(self) -> str:
        """Return the direction of the shortest strand.

        Can be 'fwd', 'rev', or 'both' if both strands are equally short.
        """
        if self.length_fwd > self.length_rev:
            return "rev"
        elif self.length_fwd < self.length_rev:
            return "fwd"
        else:
            return "both"

    @property
    def length_doublestrand(self) -> int:
        """Length of the double strand portion between the two mapping strands.

        ??? info "Click to expand for visual representations of each type"

            ```
            Ref 5' -> 3'
            Fwd 5' -> 3'
            Rev 3' <- 5'

            blunt_double (12)
                NNNNNNNNNNNN
                NNNNNNNNNNNN

            blunt_fiveprime_fwd_overhang (9)
                NNNNNNNNNNNN
                NNNNNNNNN___

            blunt_fiveprime_rev_overhang (9)
                NNNNNNNNN___
                NNNNNNNNNNNN

            blunt_threeprime_fwd_overhang (9)
                NNNNNNNNNNNN
                ___NNNNNNNNN

            blunt_threeprime_rev_overhang (9)
                ___NNNNNNNNN
                NNNNNNNNNNNN

            jagged_shifted_fiveprime (6)
                NNNNNNNNN___
                ___NNNNNNNNN

            jagged_shifted_threeprime (6)
                ___NNNNNNNNN
                NNNNNNNNN___

            jagged_encompassing_fwd_overhang (6)
                NNNNNNNNNNNN
                ___NNNNNN___

            jagged_encompassing_rev_overhang (6)
                ___NNNNNN___
                NNNNNNNNNNNN

            ```

        """
        start = max(self.fwd.pysam.reference_start, self.rev.pysam.reference_start)
        end = min(self.fwd.pysam.reference_end, self.rev.pysam.reference_end)  # ty: ignore
        return end - start

    @property
    def length_span(self) -> int:
        """Length of from the most 5' aligned based to the 3' most aligned based on
        either strand.

        ??? info "Click to expand for visual representations of each type"

            ```
            Ref 5' -> 3'
            Fwd 5' -> 3'
            Rev 3' <- 5'

            blunt_double (12)
                NNNNNNNNNNNN
                NNNNNNNNNNNN

            blunt_fiveprime_fwd_overhang (12)
                NNNNNNNNNNNN
                NNNNNNNNN___

            blunt_fiveprime_rev_overhang (12)
                NNNNNNNNN___
                NNNNNNNNNNNN

            blunt_threeprime_fwd_overhang (12)
                NNNNNNNNNNNN
                ___NNNNNNNNN

            blunt_threeprime_rev_overhang (12)
                ___NNNNNNNNN
                NNNNNNNNNNNN

            jagged_shifted_fiveprime (12)
                NNNNNNNNN___
                ___NNNNNNNNN

            jagged_shifted_threeprime (12)
                ___NNNNNNNNN
                NNNNNNNNN___

            jagged_encompassing_fwd_overhang (12)
                NNNNNNNNNNNN
                ___NNNNNN___

            jagged_encompassing_rev_overhang (12)
                ___NNNNNN___
                NNNNNNNNNNNN

            ```
        """
        start = min(self.fwd.pysam.reference_start, self.rev.pysam.reference_start)
        end = max(self.fwd.pysam.reference_end, self.rev.pysam.reference_end)  # ty: ignore
        return end - start

    @property
    def length_reference_fiveprime_fwd(self) -> int:
        """Length of the overhang of the 5' reference overhang on the forward strand"""
        return max(self.rev.pysam.reference_start - self.fwd.pysam.reference_start, 0)

    @property
    def length_reference_fiveprime_rev(self) -> int:
        """Length of the overhang of the 5' reference overhang on the reverse strand"""
        return max(self.fwd.pysam.reference_start - self.rev.pysam.reference_start, 0)

    @property
    def length_reference_threeprime_fwd(self) -> int:
        """Length of the overhang of the 3' reference overhang on the forward strand"""
        return max(self.fwd.pysam.reference_end - self.rev.pysam.reference_end, 0)  # ty: ignore

    @property
    def length_reference_threeprime_rev(self) -> int:
        """Length of the overhang of the 3' reference overhang on the reverse strand"""
        return max(self.rev.pysam.reference_end - self.fwd.pysam.reference_end, 0)  # ty: ignore

    @property
    def length_reference_fiveprime(self) -> int:
        """Length of the overhang of the 5' reference overhang, regardless of the strand"""
        return max(
            self.length_reference_fiveprime_fwd, self.length_reference_fiveprime_rev
        )

    @property
    def length_reference_threeprime(self) -> int:
        """Length of the overhang of the 3' reference overhang, regardless of the strand"""
        return max(
            self.length_reference_threeprime_fwd, self.length_reference_threeprime_rev
        )

    @property
    def length_endrepair(self) -> int:
        """Length of the molecule if end repair was done

        End repairs blunts a double stranded DNA molecule by filling the 5' overhang
        and excising single strand bases from the 3' ends.

        ??? info "Click to expand for visual representations of each type"

            ```
            Examples below show the original double strand DNA molecule before (left) and
            after (right) *in silico* end repair. Numbers between parenthesis indicate the
            calculated length (after end repair).

            Fwd 5' -> 3'
            Rev 3' <- 5'

            blunt_double (12)
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN

            blunt_fiveprime_fwd_overhang (9)
                NNNNNNNNNNNN    ->    NNNNNNNNN___
                NNNNNNNNN___    ->    NNNNNNNNN___

            blunt_fiveprime_rev_overhang (12)
                NNNNNNNNN___    ->    NNNNNNNNNNNN
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN

            blunt_threeprime_fwd_overhang (9)
                ___NNNNNNNNN    ->    ___NNNNNNNNN
                NNNNNNNNNNNN    ->    ___NNNNNNNNN

            blunt_threeprime_rev_overhang (12)
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN
                ___NNNNNNNNN    ->    NNNNNNNNNNNN

            jagged_shifted_fiveprime (12)
                NNNNNNNNN___    ->    NNNNNNNNNNNN
                ___NNNNNNNNN    ->    NNNNNNNNNNNN

            jagged_shifted_threeprime (6)
                ___NNNNNNNNN    ->    ___NNNNNN___
                NNNNNNNNN___    ->    ___NNNNNN___

            jagged_encompassing_fwd_overhang (9)
                NNNNNNNNNNNN    ->    NNNNNNNNN___
                ___NNNNNN___    ->    NNNNNNNNN___

            jagged_encompassing_rev_overhang (9)
                ___NNNNNN___    ->    ___NNNNNNNNN
                NNNNNNNNNNNN    ->    ___NNNNNNNNN
            ```
        """
        start = self.fwd.pysam.reference_start
        end = self.rev.pysam.reference_end
        return end - start  # ty: ignore

    @property
    def endrepair_removed_bases(self) -> int:
        """Number of bases that would be removed at the 3' overhangs by end repair.

        ??? info "Click to expand for visual representations of each type"

            ```
            Examples show the original double strand DNA molecule before (left) and after
            in silico end repair (right). F characters indicate filled bases. R symbols
            indicate removed bases. Numbers between parenthesis indicate the number of
            removed bases (after end repair).

            Fwd 5' -> 3'
            Rev 3' <- 5'

            blunt_double (0)
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN

            blunt_fiveprime_fwd_overhang (3)
                NNNNNNNNNNNN    ->    NNNNNNNNNRRR
                NNNNNNNNN___    ->    NNNNNNNNN___

            blunt_fiveprime_rev_overhang (0)
                NNNNNNNNN___    ->    NNNNNNNNNFFF
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN

            blunt_threeprime_fwd_overhang (3)
                ___NNNNNNNNN    ->    ___NNNNNNNNN
                NNNNNNNNNNNN    ->    RRRNNNNNNNNN

            blunt_threeprime_rev_overhang (0)
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN
                ___NNNNNNNNN    ->    FFFNNNNNNNNN

            jagged_shifted_fiveprime (0)
                NNNNNNNNN___    ->    NNNNNNNNNFFF
                ___NNNNNNNNN    ->    FFFNNNNNNNNN

            jagged_shifted_threeprime (6)
                ___NNNNNNNNN    ->    ___NNNNNNRRR
                NNNNNNNNN___    ->    RRRNNNNNN___

            jagged_encompassing_fwd_overhang (3)
                NNNNNNNNNNNN    ->    NNNNNNNNNRRR
                ___NNNNNN___    ->    FFFNNNNNN___

            jagged_encompassing_rev_overhang (3)
                ___NNNNNN___    ->    ___NNNNNNFFF
                NNNNNNNNNNNN    ->    RRRNNNNNNNNN
            ```
        """

        return (
            self.length_reference_fiveprime_rev + self.length_reference_threeprime_fwd
        )

    @property
    def endrepair_filled_bases(self) -> int:
        """Number of bases that would be filled at the 5' overhangs by end repair.

        ??? info "Click to expand for visual representations of each type"

            ```
            Examples show the original double strand DNA molecule before (left) and after
            in silico end repair (right). F characters indicate filled bases. R symbols
            indicate removed bases. Numbers between parenthesis indicate the number of
            filled bases (after end repair).

            Fwd 5' -> 3'
            Rev 3' <- 5'

            blunt_double (0)
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN

            blunt_fiveprime_fwd_overhang (0)
                NNNNNNNNNNNN    ->    NNNNNNNNNRRR
                NNNNNNNNN___    ->    NNNNNNNNN___

            blunt_fiveprime_rev_overhang (3)
                NNNNNNNNN___    ->    NNNNNNNNNFFF
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN

            blunt_threeprime_fwd_overhang (0)
                ___NNNNNNNNN    ->    ___NNNNNNNNN
                NNNNNNNNNNNN    ->    RRRNNNNNNNNN

            blunt_threeprime_rev_overhang (3)
                NNNNNNNNNNNN    ->    NNNNNNNNNNNN
                ___NNNNNNNNN    ->    FFFNNNNNNNNN

            jagged_shifted_fiveprime (6)
                NNNNNNNNN___    ->    NNNNNNNNNFFF
                ___NNNNNNNNN    ->    FFFNNNNNNNNN

            jagged_shifted_threeprime (0)
                ___NNNNNNNNN    ->    ___NNNNNNRRR
                NNNNNNNNN___    ->    RRRNNNNNN___

            jagged_encompassing_fwd_overhang (3)
                NNNNNNNNNNNN    ->    NNNNNNNNNRRR
                ___NNNNNN___    ->    FFFNNNNNN___

            jagged_encompassing_rev_overhang (3)
                ___NNNNNN___    ->    ___NNNNNNFFF
                NNNNNNNNNNNN    ->    RRRNNNNNNNNN
            ```
        """

        return (
            self.length_reference_fiveprime_fwd + self.length_reference_threeprime_rev
        )

    @property
    def name(self) -> str | None:
        """Name of the read"""
        return self.fwd._pysam.query_name

    def print_structure(self):
        """Print the dsDNA sequences, for overhang visualization"""
        fwd_str_print, rev_str_print = self.get_structure_strings()

        print(f"Fwd: {fwd_str_print}")
        print(f"Rev: {rev_str_print}")

    def get_structure_strings(self) -> Tuple[str, str]:
        """Get the forward and reverse strings, padded with '_' in the overhangs

        Returns:
            A tuple with the forward padded sequence as the first item, and the
                reverse padded sequence as the second item.
        """

        fwd_str_print = (
            "_" * self.length_reference_fiveprime_rev  # ty: ignore
            + self.fwd.pysam.query_sequence
            + "_" * self.length_reference_threeprime_rev
        )
        rev_str_print = (
            "_" * self.length_reference_fiveprime_fwd  # ty: ignore
            + self.rev.pysam.query_sequence
            + "_" * self.length_reference_threeprime_fwd
        )

        return fwd_str_print, rev_str_print

    @property
    def is_duplicate(self) -> bool:
        """Return true if read is marked as duplicate"""
        return self.fwd.pysam.is_duplicate


def to_dict(obj):
    return {
        name: getattr(obj, name)
        for name in dir(obj)
        if not name.startswith("_") and not callable(getattr(obj, name))
    }
