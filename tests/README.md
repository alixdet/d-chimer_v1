# d-chimer Test Suite

This directory contains unit tests for the d-chimer package.

## Running Tests

### Prerequisites

Install test dependencies:
```bash
pip install pytest pytest-cov
# or use conda
conda env update -f environment.yml
```

### Run All Tests

```bash
pytest
```

### Run Tests with Coverage Report

```bash
pytest --cov=dchimer --cov-report=html tests/
```

This generates an HTML coverage report in `htmlcov/index.html`.

### Run Specific Test File

```bash
pytest tests/test_input_validator.py
```

### Run Specific Test Class or Function

```bash
pytest tests/test_input_validator.py::TestValidateFastaFile::test_valid_fasta_file
```

### Run Tests with Verbose Output

```bash
pytest -v
```

### Run Tests Matching a Pattern

```bash
pytest -k "validate"  # Runs all tests with "validate" in the name
```

## Test Structure

- **test_input_validator.py**: Tests for FASTA file validation
- **test_config_validator.py**: Tests for configuration validation
- **test_classes.py**: Tests for core classes (Subject, Contig, CsvIO)
- **conftest.py**: Pytest fixtures for common test data

## Test Data

Fixtures in `conftest.py` provide:
- Sample FASTA files
- Sample BLAST CSV output (BLASTn and BLASTx)
- Sample configuration files
- Temporary directories for file operations

## Adding New Tests

1. Create test files with `test_*.py` naming
2. Write test functions/classes with `test_*` naming
3. Use fixtures from `conftest.py` as needed
4. Add pytest markers (`@pytest.mark.slow`, etc.) for organization

## Continuous Integration

For CI/CD pipelines, use:
```bash
pytest --cov=dchimer --cov-report=xml tests/
```

This generates an XML coverage report compatible with services like Codecov.
