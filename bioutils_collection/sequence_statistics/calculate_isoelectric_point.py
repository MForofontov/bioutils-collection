"""Calculate isoelectric point (pI) of proteins."""


def calculate_isoelectric_point(seq: str) -> float:
    """
    Calculate theoretical isoelectric point (pI) of a protein sequence.

    Parameters
    ----------
    seq : str
        Protein sequence (single-letter amino acid codes).

    Returns
    -------
    float
        Isoelectric point (pH value).

    Raises
    ------
    TypeError
        If seq is not a string.
    ValueError
        If seq contains invalid amino acid codes or is empty.

    Examples
    --------
    >>> calculate_isoelectric_point('ACDEFGH')
    4.06
    >>> calculate_isoelectric_point('KKKKK')
    11.15

    Notes
    -----
    Solves Henderson-Hasselbalch net charge by bisection on pH [0, 14].
    Uses fixed pKa values: N-terminus 9.69, C-terminus 2.34, C 8.33,
    D 3.86, E 4.25, H 6.00, K 10.53, R 12.48, Y 10.07. This is a
    fixed-pKa model; terminal sequence context and modifications are ignored.
    Results need not match models using other pKa sets.

    Complexity
    ----------
    Time: O(n), Space: O(1)
    """
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")

    seq = seq.upper()

    valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    if not seq:
        raise ValueError("Sequence cannot be empty")

    if not all(aa in valid_amino_acids for aa in seq):
        raise ValueError("Sequence contains invalid amino acid codes")

    # pKa values for ionizable groups
    pka_nterm = 9.69
    pka_cterm = 2.34

    pka_side_chains = {
        "C": 8.33,  # Cysteine
        "D": 3.86,  # Aspartic acid
        "E": 4.25,  # Glutamic acid
        "H": 6.00,  # Histidine
        "K": 10.53,  # Lysine
        "R": 12.48,  # Arginine
        "Y": 10.07,  # Tyrosine
    }

    counts = {aa: seq.count(aa) for aa in pka_side_chains}
    low, high = 0.0, 14.0
    # Net charge is monotone decreasing in pH. Bisect well beyond 0.01-pH
    # output precision, including terminal, cysteine and tyrosine charges.
    for _ in range(60):
        ph = (low + high) / 2
        positive = 1 / (1 + 10 ** (ph - pka_nterm))
        negative = 1 / (1 + 10 ** (pka_cterm - ph))
        for aa in ("K", "R", "H"):
            positive += counts[aa] / (1 + 10 ** (ph - pka_side_chains[aa]))
        for aa in ("D", "E", "C", "Y"):
            negative += counts[aa] / (1 + 10 ** (pka_side_chains[aa] - ph))
        if positive > negative:
            low = ph
        else:
            high = ph
    return round((low + high) / 2, 2)


__all__ = ["calculate_isoelectric_point"]
