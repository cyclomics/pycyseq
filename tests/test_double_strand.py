"""
Unit tests for DoubleStrandAlignment.

Strategy
--------
Neither Alignment nor DoubleStrandAlignment care what concrete object sits
behind `Alignment.pysam` -- they only ever read a handful of attributes off
it: reference_start, reference_end, is_forward, query_sequence, query_name,
is_mapped. `Alignment.__init__` doesn't type-check its `pysam` argument at
all, and `DoubleStrandAlignment.__init__` only checks `isinstance(x,
Alignment)`. So instead of constructing a real (annoying-to-build)
pysam.AlignedSegment, we use a small fake dataclass with just those
attributes, and wrap it in a real Alignment.
"""

from dataclasses import dataclass
from typing import Optional

import pytest

from pycyseq.alignment import Alignment
from pycyseq.double_strand import DoubleStrandAlignment, to_dict


@dataclass
class FakeAlignedSegment:
    """Minimal stand-in for pysam.AlignedSegment."""

    reference_start: int
    reference_end: int
    is_forward: bool
    query_sequence: Optional[str] = None
    query_name: str = "read1"
    is_mapped: bool = True
    is_duplicate: bool = False

    def __post_init__(self):
        if self.query_sequence is None:
            self.query_sequence = "N" * (self.reference_end - self.reference_start)


def make_alignment(start: int, end: int, is_forward: bool, **kwargs) -> Alignment:
    """Build a real Alignment wrapping a FakeAlignedSegment."""
    return Alignment(FakeAlignedSegment(start, end, is_forward, **kwargs))  # ty: ignore


def make_ds(fwd_coords, rev_coords, **kwargs) -> DoubleStrandAlignment:
    """Build a DoubleStrandAlignment from (start, end) tuples."""
    fwd = make_alignment(*fwd_coords, is_forward=True, **kwargs)
    rev = make_alignment(*rev_coords, is_forward=False, **kwargs)
    return DoubleStrandAlignment(fwd, rev)


# id, fwd (start, end), rev (start, end) -- see module docstring for diagrams
CASE_COORDS = [
    # NNN NNN NNN NNN
    # NNN NNN NNN NNN
    ("blunt_double", (0, 12), (0, 12)),
    # NNN NNN NNN NNN
    # NNN NNN NNN
    ("blunt_fiveprime_fwd_overhang", (0, 12), (0, 9)),
    # NNN NNN NNN
    # NNN NNN NNN NNN
    ("blunt_fiveprime_rev_overhang", (0, 9), (0, 12)),
    # NNN NNN NNN NNN
    #     NNN NNN NNN
    ("blunt_threeprime_fwd_overhang", (0, 12), (3, 12)),
    #     NNN NNN NNN
    # NNN NNN NNN NNN
    ("blunt_threeprime_rev_overhang", (3, 12), (0, 12)),
    # NNN NNN NNN
    #     NNN NNN NNN
    ("jagged_shifted_fiveprime", (0, 9), (3, 12)),
    #     NNN NNN NNN
    # NNN NNN NNN
    ("jagged_shifted_threeprime", (3, 12), (0, 9)),
    # NNN NNN NNN NNN
    #     NNN NNN
    ("jagged_encompassing_fwd_overhang", (0, 12), (3, 9)),
    #     NNN NNN
    # NNN NNN NNN NNN
    ("jagged_encompassing_rev_overhang", (3, 9), (0, 12)),
]


def cases_with_expected(expected):
    """Pair the canonical coordinate cases with a per-property list of
    expected values (same order as CASE_COORDS) -> list of pytest.param.

    Usage:
        @pytest.mark.parametrize(
            "fwd_coords, rev_coords, expected",
            cases_with_expected([
                None,  # blunt_double
                None,  # blunt_fiveprime_fwd_overhang
                ...
            ]),
        )
    """
    assert len(expected) == len(CASE_COORDS), (
        f"expected {len(CASE_COORDS)} values (one per case), got {len(expected)}"
    )
    return [
        pytest.param(fwd, rev, exp, id=name)
        for (name, fwd, rev), exp in zip(CASE_COORDS, expected)
    ]


# ---------------------------------------------------------------------------
# Construction / validation
# ---------------------------------------------------------------------------


def test_rejects_non_alignment_fwd():
    rev = make_alignment(0, 12, is_forward=False)
    with pytest.raises(ValueError):
        DoubleStrandAlignment(fwd="not an alignment", rev=rev)  # type: ignore[arg-type]


def test_rejects_non_alignment_rev():
    fwd = make_alignment(0, 12, is_forward=True)
    with pytest.raises(ValueError):
        DoubleStrandAlignment(fwd=fwd, rev="not an alignment")  # type: ignore[arg-type]


def test_rejects_reverse_aligned_fwd():
    fwd = make_alignment(0, 12, is_forward=False)  # wrong direction
    rev = make_alignment(0, 12, is_forward=False)
    with pytest.raises(ValueError):
        DoubleStrandAlignment(fwd, rev)


def test_rejects_forward_aligned_rev():
    fwd = make_alignment(0, 12, is_forward=True)
    rev = make_alignment(0, 12, is_forward=True)  # wrong direction
    with pytest.raises(ValueError):
        DoubleStrandAlignment(fwd, rev)


def test_fwd_rev_accessors_return_the_given_alignments():
    fwd = make_alignment(0, 12, is_forward=True)
    rev = make_alignment(0, 12, is_forward=False)
    ds = DoubleStrandAlignment(fwd, rev)
    assert ds.fwd is fwd
    assert ds.rev is rev


# ---------------------------------------------------------------------------
# `type` classification
#
# Coordinates below come straight from the ASCII diagrams in the docstrings
# (each block is 12 bases wide, overhangs are 3 or 6 bases).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected_type",
    [
        ((0, 12), (0, 12), "blunt_double"),
        ((0, 12), (0, 9), "blunt_fiveprime_fwd_overhang"),
        ((0, 9), (0, 12), "blunt_fiveprime_rev_overhang"),
        ((0, 12), (3, 12), "blunt_threeprime_fwd_overhang"),
        ((3, 12), (0, 12), "blunt_threeprime_rev_overhang"),
        ((0, 9), (3, 12), "jagged_shifted_fiveprime"),
        ((3, 12), (0, 9), "jagged_shifted_threeprime"),
        ((0, 12), (3, 9), "jagged_encompassing_fwd_overhang"),
        ((3, 9), (0, 12), "jagged_encompassing_rev_overhang"),
    ],
    ids=[name for name, _, _ in CASE_COORDS],
)
def test_type_classification(fwd_coords, rev_coords, expected_type):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.type == expected_type


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            12,  # blunt_fiveprime_fwd_overhang
            9,  # blunt_fiveprime_rev_overhang
            12,  # blunt_threeprime_fwd_overhang
            9,  # blunt_threeprime_rev_overhang
            9,  # jagged_shifted_fiveprime
            9,  # jagged_shifted_threeprime
            12,  # jagged_encompassing_fwd_overhang
            6,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_fwd(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_fwd == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            9,  # blunt_fiveprime_fwd_overhang
            12,  # blunt_fiveprime_rev_overhang
            9,  # blunt_threeprime_fwd_overhang
            12,  # blunt_threeprime_rev_overhang
            9,  # jagged_shifted_fiveprime
            9,  # jagged_shifted_threeprime
            6,  # jagged_encompassing_fwd_overhang
            12,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_rev(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_rev == expected


# # ---------------------------------------------------------------------------
# # length_sum / length_long / length_short / longest_strand / shortest_strand
# # ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            24,  # blunt_double
            21,  # blunt_fiveprime_fwd_overhang
            21,  # blunt_fiveprime_rev_overhang
            21,  # blunt_threeprime_fwd_overhang
            21,  # blunt_threeprime_rev_overhang
            18,  # jagged_shifted_fiveprime
            18,  # jagged_shifted_threeprime
            18,  # jagged_encompassing_fwd_overhang
            18,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_sum(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_sum == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            12,  # blunt_fiveprime_fwd_overhang
            12,  # blunt_fiveprime_rev_overhang
            12,  # blunt_threeprime_fwd_overhang
            12,  # blunt_threeprime_rev_overhang
            9,  # jagged_shifted_fiveprime
            9,  # jagged_shifted_threeprime
            12,  # jagged_encompassing_fwd_overhang
            12,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_long(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_long == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            9,  # blunt_fiveprime_fwd_overhang
            9,  # blunt_fiveprime_rev_overhang
            9,  # blunt_threeprime_fwd_overhang
            9,  # blunt_threeprime_rev_overhang
            9,  # jagged_shifted_fiveprime
            9,  # jagged_shifted_threeprime
            6,  # jagged_encompassing_fwd_overhang
            6,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_short(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_short == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            "both",  # blunt_double
            "fwd",  # blunt_fiveprime_fwd_overhang
            "rev",  # blunt_fiveprime_rev_overhang
            "fwd",  # blunt_threeprime_fwd_overhang
            "rev",  # blunt_threeprime_rev_overhang
            "both",  # jagged_shifted_fiveprime
            "both",  # jagged_shifted_threeprime
            "fwd",  # jagged_encompassing_fwd_overhang
            "rev",  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_longest_strand(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.longest_strand == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            "both",  # blunt_double
            "rev",  # blunt_fiveprime_fwd_overhang
            "fwd",  # blunt_fiveprime_rev_overhang
            "rev",  # blunt_threeprime_fwd_overhang
            "fwd",  # blunt_threeprime_rev_overhang
            "both",  # jagged_shifted_fiveprime
            "both",  # jagged_shifted_threeprime
            "rev",  # jagged_encompassing_fwd_overhang
            "fwd",  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_shortest_strand(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.shortest_strand == expected


# # ---------------------------------------------------------------------------
# # length_doublestrand / length_span
# # ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            9,  # blunt_fiveprime_fwd_overhang
            9,  # blunt_fiveprime_rev_overhang
            9,  # blunt_threeprime_fwd_overhang
            9,  # blunt_threeprime_rev_overhang
            6,  # jagged_shifted_fiveprime
            6,  # jagged_shifted_threeprime
            6,  # jagged_encompassing_fwd_overhang
            6,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_doublestrand(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_doublestrand == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            12,  # blunt_fiveprime_fwd_overhang
            12,  # blunt_fiveprime_rev_overhang
            12,  # blunt_threeprime_fwd_overhang
            12,  # blunt_threeprime_rev_overhang
            12,  # jagged_shifted_fiveprime
            12,  # jagged_shifted_threeprime
            12,  # jagged_encompassing_fwd_overhang
            12,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_span(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_span == expected


# ---------------------------------------------------------------------------
# length_reference_fiveprime_fwd / _rev, length_reference_threeprime_fwd / _rev
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            0,  # blunt_fiveprime_fwd_overhang
            0,  # blunt_fiveprime_rev_overhang
            3,  # blunt_threeprime_fwd_overhang
            0,  # blunt_threeprime_rev_overhang
            3,  # jagged_shifted_fiveprime
            0,  # jagged_shifted_threeprime
            3,  # jagged_encompassing_fwd_overhang
            0,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_reference_fiveprime_fwd(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_reference_fiveprime_fwd == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            0,  # blunt_fiveprime_fwd_overhang
            0,  # blunt_fiveprime_rev_overhang
            0,  # blunt_threeprime_fwd_overhang
            3,  # blunt_threeprime_rev_overhang
            0,  # jagged_shifted_fiveprime
            3,  # jagged_shifted_threeprime
            0,  # jagged_encompassing_fwd_overhang
            3,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_reference_fiveprime_rev(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_reference_fiveprime_rev == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            3,  # blunt_fiveprime_fwd_overhang
            0,  # blunt_fiveprime_rev_overhang
            0,  # blunt_threeprime_fwd_overhang
            0,  # blunt_threeprime_rev_overhang
            0,  # jagged_shifted_fiveprime
            3,  # jagged_shifted_threeprime
            3,  # jagged_encompassing_fwd_overhang
            0,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_reference_threeprime_fwd(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_reference_threeprime_fwd == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            0,  # blunt_fiveprime_fwd_overhang
            3,  # blunt_fiveprime_rev_overhang
            0,  # blunt_threeprime_fwd_overhang
            0,  # blunt_threeprime_rev_overhang
            3,  # jagged_shifted_fiveprime
            0,  # jagged_shifted_threeprime
            0,  # jagged_encompassing_fwd_overhang
            3,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_reference_threeprime_rev(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_reference_threeprime_rev == expected


# ---------------------------------------------------------------------------
# length_reference_fiveprime / length_reference_threeprime (strand-agnostic)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            0,  # blunt_fiveprime_fwd_overhang
            0,  # blunt_fiveprime_rev_overhang
            3,  # blunt_threeprime_fwd_overhang
            3,  # blunt_threeprime_rev_overhang
            3,  # jagged_shifted_fiveprime
            3,  # jagged_shifted_threeprime
            3,  # jagged_encompassing_fwd_overhang
            3,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_reference_fiveprime(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_reference_fiveprime == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            3,  # blunt_fiveprime_fwd_overhang
            3,  # blunt_fiveprime_rev_overhang
            0,  # blunt_threeprime_fwd_overhang
            0,  # blunt_threeprime_rev_overhang
            3,  # jagged_shifted_fiveprime
            3,  # jagged_shifted_threeprime
            3,  # jagged_encompassing_fwd_overhang
            3,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_reference_threeprime(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_reference_threeprime == expected


# ---------------------------------------------------------------------------
# length_endrepair / endrepair_removed_bases / endrepair_filled_bases
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            12,  # blunt_double
            9,  # blunt_fiveprime_fwd_overhang
            12,  # blunt_fiveprime_rev_overhang
            12,  # blunt_threeprime_fwd_overhang
            9,  # blunt_threeprime_rev_overhang
            12,  # jagged_shifted_fiveprime
            6,  # jagged_shifted_threeprime
            9,  # jagged_encompassing_fwd_overhang
            9,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_length_endrepair(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.length_endrepair == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            3,  # blunt_fiveprime_fwd_overhang
            0,  # blunt_fiveprime_rev_overhang
            0,  # blunt_threeprime_fwd_overhang
            3,  # blunt_threeprime_rev_overhang
            0,  # jagged_shifted_fiveprime
            6,  # jagged_shifted_threeprime
            3,  # jagged_encompassing_fwd_overhang
            3,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_endrepair_removed_bases(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.endrepair_removed_bases == expected


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            0,  # blunt_double
            0,  # blunt_fiveprime_fwd_overhang
            3,  # blunt_fiveprime_rev_overhang
            3,  # blunt_threeprime_fwd_overhang
            0,  # blunt_threeprime_rev_overhang
            6,  # jagged_shifted_fiveprime
            0,  # jagged_shifted_threeprime
            3,  # jagged_encompassing_fwd_overhang
            3,  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_endrepair_filled_bases(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.endrepair_filled_bases == expected


# ---------------------------------------------------------------------------
# get_structure_strings -> (fwd_str, rev_str)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fwd_coords, rev_coords, expected",
    cases_with_expected(
        [
            (
                "NNNNNNNNNNNN",
                "NNNNNNNNNNNN",
            ),  # blunt_double            -> (fwd_str, rev_str)
            ("NNNNNNNNNNNN", "NNNNNNNNN___"),  # blunt_fiveprime_fwd_overhang
            ("NNNNNNNNN___", "NNNNNNNNNNNN"),  # blunt_fiveprime_rev_overhang
            ("NNNNNNNNNNNN", "___NNNNNNNNN"),  # blunt_threeprime_fwd_overhang
            ("___NNNNNNNNN", "NNNNNNNNNNNN"),  # blunt_threeprime_rev_overhang
            ("NNNNNNNNN___", "___NNNNNNNNN"),  # jagged_shifted_fiveprime
            ("___NNNNNNNNN", "NNNNNNNNN___"),  # jagged_shifted_threeprime
            ("NNNNNNNNNNNN", "___NNNNNN___"),  # jagged_encompassing_fwd_overhang
            ("___NNNNNN___", "NNNNNNNNNNNN"),  # jagged_encompassing_rev_overhang
        ]
    ),
)
def test_get_structure_strings(fwd_coords, rev_coords, expected):
    ds = make_ds(fwd_coords, rev_coords)
    assert ds.get_structure_strings() == expected


# ---------------------------------------------------------------------------
# name
# ---------------------------------------------------------------------------


def test_name_comes_from_fwd_alignment():
    ds = make_ds((0, 12), (0, 12), query_name="my-read")
    assert ds.name == "my-read"


# ---------------------------------------------------------------------------
# __str__  (currently broken -- see module docstring, issue 2)
# ---------------------------------------------------------------------------


def test_str_raises_due_to_missing_aln_attribute():
    """DoubleStrandAlignment.__str__ calls self._fwd.aln, but Alignment only
    exposes `.pysam`. This documents the current (broken) behavior -- once
    `__str__` is fixed to use `.pysam` (or something sensible), replace this
    with an assertion on the actual expected string.
    """
    ds = make_ds((0, 12), (0, 12))
    with pytest.raises(AttributeError):
        str(ds)


# ---------------------------------------------------------------------------
# to_dict (module-level helper)
# ---------------------------------------------------------------------------


def test_to_dict_includes_all_public_properties():
    ds = make_ds((0, 9), (3, 12), query_name="my-read")
    d = to_dict(ds)
    assert d["name"] == "my-read"
    assert d["length_fwd"] == ds.length_fwd
    assert d["length_rev"] == ds.length_rev
    assert d["type"] == ds.type
