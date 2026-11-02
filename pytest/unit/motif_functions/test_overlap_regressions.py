from bioutils_collection import sequence_pattern_match


def test_iupac_matches_include_overlaps() -> None:
    assert sequence_pattern_match("AAAA", "AA") == [0, 1, 2]
    assert sequence_pattern_match("AAAA", "AN") == [0, 1, 2]
    assert sequence_pattern_match("AAAA", "AA", use_iupac=False) == [0, 1, 2]
