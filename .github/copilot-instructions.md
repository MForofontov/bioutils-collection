# GitHub Copilot Instructions for Bioutils Collection

## Project Overview
Bioutils Collection is a production-ready bioinformatics toolkit with 77+ optimized functions for sequence analysis, alignment, annotation, and molecular biology workflows. This codebase follows strict typing and testing standards.

## Code Style & Standards

### Type Annotations
- **REQUIRED**: All functions must have complete type annotations
- Use Python 3.10+ syntax: `dict[str, str]`, `list[int]` (not `Dict`, `List`)
- Use `|` for unions instead of `Union`: `str | int` not `Union[str, int]`
- Avoid `typing.Optional`, use `Type | None` instead
- Return types must be explicit, never use `Any` in return types
- Use `from collections.abc import Iterator` for generator types

### Function Structure
Every function must follow this exact structure:

```python
"""Brief one-line description.

Extended description explaining what the function does,
its purpose, and any important algorithmic details.

Parameters
----------
param_name : type
    Description of parameter.
param2 : type, optional
    Description with (default: value).

Returns
-------
return_type
    Description of what is returned.

Raises
------
ErrorType
    When this error occurs.

Examples
--------
>>> function_name(example_input)
expected_output
>>> function_name(another_example)
another_output

Notes
-----
Additional information, references, algorithms used.

Complexity
----------
Time: O(notation), Space: O(notation) with explanation
"""
```

### Input Validation Pattern
**ALWAYS** validate inputs in this exact order:

```python
# 1. Type validation
if not isinstance(param, expected_type):
    raise TypeError(f"param must be {expected_type.__name__}, got {type(param).__name__}")

# 2. Empty/None checks
if not param:
    return default_or_raise

# 3. Value range validation
if param < min_value or param > max_value:
    raise ValueError(f"param must be between {min_value} and {max_value}")

# 4. Content validation
if invalid_condition:
    raise ValueError("Specific error message")
```

### Module Organization
Each module follows this pattern:

```python
"""Module description."""

# Standard library imports
import math
from collections.abc import Iterator

# Third-party imports (if needed)
import numpy as np

# Local imports (avoid circular dependencies)
from .submodule import function


def public_function():
    """Public API function."""
    pass


def _private_function():
    """Internal helper function."""
    pass


__all__ = ["public_function"]
```

## Project Structure

### Module Categories
1. **alignment_functions**: Pairwise/multiple sequence alignment (Needleman-Wunsch, Smith-Waterman)
2. **annotation_functions**: BED/GFF/GTF/VCF parsing and conversion
3. **clustering_functions**: Motif clustering algorithms
4. **data_validation**: Sequence validation (DNA/RNA/protein)
5. **fasta_misc**: FASTA/FASTQ file operations, parsing, conversion
6. **gc_functions**: GC content analysis and profiling
7. **motif_functions**: Pattern matching and consensus generation
8. **repeat_functions**: Tandem repeats and palindrome detection
9. **restriction_functions**: Restriction enzyme site identification
10. **sequence_operations**: Core sequence manipulation (reverse complement, ORFs, CpG islands)
11. **sequence_statistics**: Statistical analysis (codon usage, entropy, complexity)
12. **translation_functions**: DNA↔RNA↔Protein translation

### File Naming Convention
- Use snake_case for all files: `needleman_wunsch.py`
- One primary function per file (function name = file name)
- Module `__init__.py` exports all public functions

## Common Patterns

### Sequence Validation
```python
seq = seq.upper()
valid_bases = set("ATGC")
if not all(base in valid_bases for base in seq):
    raise ValueError("Sequence contains invalid DNA bases")
```

### Iterator/Generator Pattern
```python
from collections.abc import Iterator

def parse_format(data: str) -> Iterator[tuple[str, str]]:
    """Parse and yield results."""
    for item in data.splitlines():
        if condition:
            yield (header, sequence)
```

### Dictionary Return Pattern
```python
from typing import Any

def analyze_sequence(seq: str) -> dict[str, Any]:
    """Return structured analysis results."""
    return {
        "is_valid": bool_value,
        "length": int_value,
        "statistics": float_value,
    }
```

## Testing Requirements

### Test Structure and Organization

#### Directory Structure
Tests mirror the package structure exactly:
```
pytest/
├── conftest.py                    # Shared fixtures and configuration
└── unit/
    ├── alignment_functions/
    │   ├── test_needleman_wunsch.py
    │   ├── test_smith_waterman.py
    │   └── ...
    ├── annotation_functions/
    │   ├── test_parse_bed.py
    │   ├── test_parse_gff.py
    │   └── ...
    └── [other_modules]/
        └── test_[function_name].py
```

#### File Naming Convention
- **Pattern**: `test_[function_name].py`
- **Location**: `pytest/unit/[module_name]/test_[function_name].py`
- **One test file per function** (matches source file)

#### Test Function Naming
```python
# Pattern: test_[function_name]_[scenario]_[expected_behavior]
test_needleman_wunsch_identical_sequences_returns_perfect_score()
test_needleman_wunsch_empty_sequence_raises_value_error()
test_parse_fasta_valid_input_yields_correct_pairs()
test_gc_content_all_gc_returns_hundred_percent()
```

### Test Categories (Required Coverage)

Every function must have tests for:

#### 1. Normal Functionality Tests
```python
def test_function_name_basic_case() -> None:
    """Test case 1: Test function with typical valid input."""
    result = function("ATGC")
    assert result == expected_value
    assert isinstance(result, ExpectedType)

def test_function_name_with_options() -> None:
    """Test case 2: Test function with non-default parameters."""
    result = function("ATGC", option=True, threshold=0.5)
    assert result == expected_with_options

def test_function_name_multiple_inputs() -> None:
    """Test case 3: Test function with multiple valid inputs."""
    result = function("ATGC", "GGCC")
    assert result is not None
```

#### 2. Edge Case Tests
```python
def test_function_name_empty_input() -> None:
    """Test case 4: Test function with empty string."""
    result = function("")
    assert result == ""  # or appropriate default

def test_function_name_single_character() -> None:
    """Test case 5: Test function with minimal valid input."""
    result = function("A")
    assert result is not None

def test_function_name_very_long_input() -> None:
    """Test case 6: Test function with large input."""
    long_seq = "ATGC" * 10000
    result = function(long_seq)
    assert len(result) == expected_length

def test_function_name_boundary_values() -> None:
    """Test case 7: Test function at parameter boundaries."""
    result = function("ATGC", window=1)  # Minimum
    assert result is not None
    result = function("ATGC", threshold=1.0)  # Maximum
    assert result is not None

def test_function_name_lowercase() -> None:
    """Test case 8: Test function handles lowercase input."""
    result_upper = function("ATGC")
    result_lower = function("atgc")
    assert result_upper == result_lower
```

#### 3. Error Handling Tests (Critical)
```python
def test_function_name_type_error_not_string() -> None:
    """Test case 9: TypeError when input is not a string."""
    with pytest.raises(TypeError, match="must be str, got int"):
        function(12345)  # int instead of str

def test_function_name_type_error_none() -> None:
    """Test case 10: TypeError when input is None."""
    with pytest.raises(TypeError, match="must be str, got NoneType"):
        function(None)

def test_function_name_value_error_invalid_chars() -> None:
    """Test case 11: ValueError for invalid characters."""
    with pytest.raises(ValueError, match="Invalid DNA bases"):
        function("ATGCX")  # Invalid character

def test_function_name_value_error_negative() -> None:
    """Test case 12: ValueError for negative parameter."""
    with pytest.raises(ValueError, match="must be positive"):
        function("ATGC", window=-1)

def test_function_name_value_error_empty() -> None:
    """Test case 13: ValueError when empty input not allowed."""
    with pytest.raises(ValueError, match="cannot be empty"):
        function("")
```

#### 4. Data Validation Tests
```python
def test_function_name_validates_sequence_content() -> None:
    """Test case 14: Validate DNA/RNA/protein content."""
    with pytest.raises(ValueError, match="Invalid.*bases"):
        function("ATGCXYZ")  # Invalid bases

def test_function_name_handles_mixed_case() -> None:
    """Test case 15: Handle mixed case input correctly."""
    result = function("AtGc")
    assert result == expected_value

def test_function_name_preserves_case() -> None:
    """Test case 16: Verify case preservation when expected."""
    result = function("AtGc")
    assert "A" in result  # Upper preserved
    assert "t" in result  # Lower preserved
```

#### 5. Return Type Tests
```python
def test_function_name_returns_correct_type() -> None:
    """Test case 17: Verify function returns documented type."""
    result = function("ATGC")
    assert isinstance(result, str)  # or expected type

def test_function_name_dict_return_has_required_keys() -> None:
    """Test case 18: Dictionary return contains all documented keys."""
    result = function("ATGC")
    assert isinstance(result, dict)
    assert "is_valid" in result
    assert "length" in result
    assert "gc_content" in result

def test_function_name_iterator_yields_correct_types() -> None:
    """Test case 19: Iterator yields documented types."""
    results = list(function("ATGC"))
    assert all(isinstance(item, tuple) for item in results)
    assert all(len(item) == 2 for item in results)
```

### Pytest Markers

Use markers to categorize tests:

```python
import pytest

@pytest.mark.unit
@pytest.mark.alignment
def test_needleman_wunsch_basic():
    """Basic alignment test."""
    pass

@pytest.mark.unit
@pytest.mark.slow
@pytest.mark.alignment
def test_needleman_wunsch_large_sequences():
    """Test with large sequences (slow)."""
    pass

@pytest.mark.unit
@pytest.mark.fasta
def test_parse_fasta_large_file():
    """Test FASTA parser with large file."""
    pass
```

**Available Markers**:
- `@pytest.mark.unit` - All unit tests (required)
- `@pytest.mark.slow` - Tests taking >1 second
- `@pytest.mark.alignment` - Alignment function tests
- `@pytest.mark.annotation` - Annotation function tests
- `@pytest.mark.fasta` - FASTA/FASTQ tests
- `@pytest.mark.translation` - Translation tests
- `@pytest.mark.validation` - Validation tests
- Module-specific markers for other categories

### Test Fixtures

#### Common Fixtures (in `conftest.py`)

```python
import pytest

@pytest.fixture
def dna_sequence():
    """Standard DNA sequence for testing."""
    return "ATGCATGCATGC"

@pytest.fixture
def protein_sequence():
    """Standard protein sequence for testing."""
    return "MHACHM"

@pytest.fixture
def fasta_data():
    """Sample FASTA formatted data."""
    return ">seq1\nATGC\n>seq2\nGGTT\n"

@pytest.fixture
def bed_data():
    """Sample BED format data."""
    return "chr1\t1000\t2000\tfeature1\t100\t+\n"

@pytest.fixture
def long_sequence():
    """Long sequence for performance testing."""
    return "ATGC" * 10000

@pytest.fixture
def invalid_sequence():
    """Sequence with invalid characters."""
    return "ATGCXYZ123"
```

#### Using Fixtures in Tests

```python
def test_function_with_fixture(dna_sequence):
    """Test using fixture data."""
    result = function(dna_sequence)
    assert result is not None

def test_function_with_multiple_fixtures(dna_sequence, long_sequence):
    """Test using multiple fixtures."""
    result1 = function(dna_sequence)
    result2 = function(long_sequence)
    assert len(result2) > len(result1)
```

### Parametrized Tests

Use parametrization for testing multiple inputs:

```python
@pytest.mark.parametrize("sequence,expected", [
    ("ATGC", "GCAT"),
    ("AAAA", "TTTT"),
    ("CCCC", "GGGG"),
    ("ATCG", "TAGC"),
])
def test_reverse_complement_various_inputs(sequence, expected):
    """Test reverse complement with various inputs."""
    assert reverse_complement(sequence) == expected

@pytest.mark.parametrize("invalid_input", [
    12345,
    None,
    [],
    {},
    12.34,
])
def test_function_rejects_non_string_types(invalid_input):
    """Test function rejects various non-string types."""
    with pytest.raises(TypeError):
        function(invalid_input)

@pytest.mark.parametrize("sequence,gc", [
    ("ATGC", 50.0),
    ("AAAA", 0.0),
    ("GGGG", 100.0),
    ("ATATATAT", 0.0),
])
def test_gc_content_exact_values(sequence, gc):
    """Test GC content calculation with known values."""
    assert gc_content(sequence) == gc
```

### Assertion Patterns

#### Value Assertions
```python
# Exact equality
assert result == expected

# Approximate equality (floats)
assert abs(result - expected) < 1e-6
assert pytest.approx(result) == expected

# Type checking
assert isinstance(result, ExpectedType)
assert type(result) is ExpectedType

# Container assertions
assert len(result) == expected_length
assert item in result
assert result[0] == first_item
```

#### Exception Assertions
```python
# Basic exception check
with pytest.raises(ValueError):
    function(invalid_input)

# Check exception message
with pytest.raises(ValueError) as exc_info:
    function(invalid_input)
assert "expected substring" in str(exc_info.value)

# Check exception attributes
with pytest.raises(ValueError) as exc_info:
    function(invalid_input)
assert exc_info.value.args[0] == "Expected message"
```

#### Iterator/Generator Assertions
```python
# Convert to list and check
results = list(function(input_data))
assert len(results) == expected_count
assert results[0] == first_expected

# Check all items
assert all(isinstance(item, tuple) for item in results)
assert all(len(item) == 2 for item in results)
```

### Test Data Management

#### Inline Test Data (Small)
```python
def test_with_inline_data():
    """Test with small inline data."""
    data = "ATGC"
    result = function(data)
    assert result == "GCAT"
```

#### Fixture Test Data (Reusable)
```python
@pytest.fixture
def test_data():
    return "ATGCATGC"

def test_with_fixture_data(test_data):
    result = function(test_data)
    assert result is not None
```

#### External Test Data (Large)
```python
# For large test files (not in repo)
@pytest.fixture
def large_fasta_file(tmp_path):
    """Create temporary large FASTA file."""
    file_path = tmp_path / "test.fasta"
    with open(file_path, 'w') as f:
        for i in range(1000):
            f.write(f">seq{i}\n")
            f.write("ATGC" * 100 + "\n")
    return str(file_path)
```

### Coverage Expectations

**Minimum Requirements**:
- **Line Coverage**: >90% for all modules
- **Branch Coverage**: >85% for all modules
- **Function Coverage**: 100% (all public functions tested)

**Measure Coverage**:
```bash
pytest --cov=bioutils_collection --cov-report=html
pytest --cov=bioutils_collection --cov-report=term-missing
```

**Coverage Exclusions** (add to source):
```python
if TYPE_CHECKING:  # pragma: no cover
    from typing import ...

if __name__ == "__main__":  # pragma: no cover
    main()
```

### Performance Testing

```python
import time
import pytest

@pytest.mark.slow
def test_function_performance_large_input():
    """Test function completes in reasonable time."""
    large_seq = "ATGC" * 100000
    
    start = time.time()
    result = function(large_seq)
    elapsed = time.time() - start
    
    assert elapsed < 5.0  # Should complete in <5 seconds
    assert result is not None

@pytest.mark.slow
def test_function_memory_efficiency():
    """Test function uses generators for large data."""
    large_data = "\n".join([f">seq{i}\nATGC" for i in range(10000)])
    
    # Should return generator/iterator, not list
    result = function(large_data)
    assert hasattr(result, '__iter__')
    assert hasattr(result, '__next__')
```

### Testing Generators/Iterators

```python
def test_generator_yields_correct_count():
    """Test generator yields expected number of items."""
    result = list(parse_fasta(fasta_data))
    assert len(result) == 2

def test_generator_yields_correct_format():
    """Test generator yields correct tuple format."""
    results = list(parse_fasta(fasta_data))
    for header, sequence in results:
        assert isinstance(header, str)
        assert isinstance(sequence, str)
        assert len(sequence) > 0

def test_generator_can_be_consumed_multiple_times():
    """Test generator behavior."""
    gen1 = parse_fasta(fasta_data)
    list1 = list(gen1)
    
    # Should create new generator
    gen2 = parse_fasta(fasta_data)
    list2 = list(gen2)
    
    assert list1 == list2
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific module
pytest pytest/unit/alignment_functions/

# Run specific test file
pytest pytest/unit/alignment_functions/test_needleman_wunsch.py

# Run specific test
pytest pytest/unit/alignment_functions/test_needleman_wunsch.py::test_needleman_wunsch_basic

# Run with markers
pytest -m unit
pytest -m "unit and not slow"
pytest -m "alignment or fasta"

# Run with coverage
pytest --cov=bioutils_collection --cov-report=html

# Run in parallel (if pytest-xdist installed)
pytest -n auto

# Verbose output
pytest -v
pytest -vv  # Very verbose

# Show print statements
pytest -s

# Stop on first failure
pytest -x

# Run last failed tests
pytest --lf
```

## Performance Considerations

### When to Use Numba
- Sequences > 30K bases: Use `@njit` decorator
- Include fallback for missing Numba
- Document performance characteristics

### Memory Efficiency
- Use generators for large file parsing
- Avoid loading entire files into memory
- Use sliding windows for large sequences

## Error Messages

### Format Standards
```python
# Type errors - show actual type
raise TypeError(f"param must be str, got {type(param).__name__}")

# Value errors - be specific
raise ValueError(f"Sequence length must be multiple of 3, got {len(seq)}")

# Invalid characters - show what's invalid
raise ValueError(f"Invalid DNA bases found: {', '.join(sorted(invalid_bases))}")
```

## Dependencies
- **Core**: numpy, scipy, scikit-learn, numba
- **Dev**: pytest, mypy, ruff
- **Python**: 3.10+ required for modern type syntax

## Known Issues to Avoid

### Type Annotation Issues
❌ **Don't**: Use old-style typing
```python
from typing import Dict, List, Optional
def func(x: Optional[Dict[str, List[int]]]) -> None:
```

✅ **Do**: Use modern syntax
```python
def func(x: dict[str, list[int]] | None) -> None:
```

### Validation Issues
❌ **Don't**: Skip type validation
```python
def func(seq: str):
    return seq.upper()  # Crashes if seq is not string
```

✅ **Do**: Always validate
```python
def func(seq: str):
    if not isinstance(seq, str):
        raise TypeError(f"seq must be str, got {type(seq).__name__}")
    return seq.upper()
```

### Return Type Issues
❌ **Don't**: Return inconsistent types
```python
def func(x: int) -> str:
    if x > 0:
        return "positive"
    return 0  # Wrong type!
```

✅ **Do**: Match declared return type
```python
def func(x: int) -> str:
    if x > 0:
        return "positive"
    return "non-positive"
```

## Complexity Documentation
Always document time and space complexity:
- `O(1)` - Constant
- `O(n)` - Linear in sequence length
- `O(n*m)` - Dynamic programming matrices
- `O(n log n)` - Sorting operations
- `O(k)` - Proportional to result size

## Special Patterns

### Genetic Code Tables
For translation functions, import from centralized module:
```python
from .genetic_code_tables import get_codon_table

codon_table = get_codon_table(table)  # Supports string, int, or dict
```

### Window/Sliding Pattern
```python
for i in range(len(seq) - window + 1):
    subseq = seq[i:i + window]
    # Process window
```

### Codon Extraction
```python
codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
```

## Ruff Configuration
Code must pass ruff checks:
- Line length: 88 characters
- Selected rules: E, F, I, UP, B, TID
- Ignored: E501 (line length - handled by formatter)

## Mypy Strict Mode
All code must pass mypy strict checks:
- No implicit `Any` types
- All functions typed
- No untyped decorators (use `# type: ignore[misc]` for Numba)

## Publishing & Versioning
- Version in `bioutils_collection/_version.py`
- SemVer: MAJOR.MINOR.PATCH
- Update CHANGELOG.md for all changes

## Quick Reference: Common Function Signatures

```python
# Simple transformation
def transform_sequence(seq: str) -> str:

# Analysis with options
def analyze_sequence(seq: str, option: int = 10) -> dict[str, Any]:

# Parser/generator
def parse_format(data: str) -> Iterator[tuple[str, str]]:

# Pairwise comparison
def compare_sequences(seq1: str, seq2: str, threshold: float = 0.5) -> float:

# File operations
def process_file(filepath: str, output: str | None = None) -> None:

# Statistical analysis
def calculate_statistics(seq: str) -> dict[str, float]:
```

## When Contributing
1. Follow the validation pattern exactly
2. Use NumPy-style docstrings with all sections
3. Add comprehensive tests (normal, edge, error cases)
4. Document complexity
5. Run: `pytest && mypy . && ruff check .`
6. Ensure type annotations are complete

---
**Remember**: This is a production library. Code quality, testing, and documentation are not optional.
