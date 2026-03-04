"""Unit tests for dchimer_config_validator module."""

import pytest
import os
from dchimer.dchimer_config_validator import (
    load_and_validate_config,
    validate_database_paths,
)


class TestLoadAndValidateConfig:
    """Tests for load_and_validate_config function."""

    def test_valid_config_file(self, sample_config_file):
        """Test loading and validating a valid config file."""
        config = load_and_validate_config(sample_config_file)
        assert config is not None
        assert 'blastn_parameters' in config
        assert 'blastx_parameters' in config
        assert config['blastn_parameters']['evalue_nt'] == 0.01

    def test_nonexistent_config_file(self):
        """Test error handling for non-existent config file."""
        with pytest.raises(FileNotFoundError):
            load_and_validate_config("/nonexistent/path/config.yaml")

    def test_missing_required_section(self, temp_dir):
        """Test error handling for missing required config section."""
        config_file = os.path.join(temp_dir, "incomplete_config.yaml")
        with open(config_file, "w") as f:
            f.write("blastn_parameters:\n  dbpath_nt: /path/to/db\n")
        
        with pytest.raises(ValueError, match="Missing configuration section"):
            load_and_validate_config(config_file)

    def test_missing_required_parameter(self, temp_dir):
        """Test error handling for missing required parameter."""
        config_file = os.path.join(temp_dir, "missing_param.yaml")
        with open(config_file, "w") as f:
            f.write("""blastn_parameters:
  dbpath_nt: /path/to/db
  nb_threads_bn: 4

filter_blastn_parameters:
  d: 10
  l: 50

blastx_parameters:
  dbpath_vrl: /path/to/viral
  dbpath_nr: /path/to/nr
  nb_threads_bx: 4
  evalue_vir: 0.1
  evalue_nr: 0.01

filter_blastx_parameters:
  d: 10
  l: 17

add_taxo_parameters:
  tax_lineages_file: /path/to/tax

blast_path: /usr/bin
""")
        
        # Missing evalue_nt in blastn_parameters
        with pytest.raises(ValueError, match="Missing parameter 'evalue_nt'"):
            load_and_validate_config(config_file)

    def test_invalid_yaml_syntax(self, temp_dir):
        """Test error handling for invalid YAML syntax."""
        config_file = os.path.join(temp_dir, "invalid_yaml.yaml")
        with open(config_file, "w") as f:
            f.write("invalid: yaml: syntax: here\n  bad indentation")
        
        with pytest.raises(ValueError, match="Invalid YAML"):
            load_and_validate_config(config_file)


class TestValidateDatabasePaths:
    """Tests for validate_database_paths function."""

    def test_existing_database_paths(self, temp_dir):
        """Test validation with existing database paths."""
        # Create mock database files
        nt_db = os.path.join(temp_dir, "nt.nin")
        tax_file = os.path.join(temp_dir, "tax.dmp")
        open(nt_db, "w").close()
        open(tax_file, "w").close()
        
        config = {
            'blastn_parameters': {'dbpath_nt': nt_db},
            'blastx_parameters': {
                'dbpath_vrl': os.path.join(temp_dir, "vir.pin"),
                'dbpath_nr': os.path.join(temp_dir, "nr.pin"),
            },
            'add_taxo_parameters': {'tax_lineages_file': tax_file},
        }
        
        # Create viral DB files
        open(os.path.join(temp_dir, "vir.pin"), "w").close()
        open(os.path.join(temp_dir, "nr.pin"), "w").close()
        
        # Should not raise if all paths exist
        validate_database_paths(config)

    def test_missing_database_files(self, temp_dir):
        """Test error handling for missing database files."""
        config = {
            'blastn_parameters': {
                'dbpath_nt': os.path.join(temp_dir, "nonexistent_nt")
            },
            'blastx_parameters': {
                'dbpath_vrl': os.path.join(temp_dir, "nonexistent_vir"),
                'dbpath_nr': os.path.join(temp_dir, "nonexistent_nr"),
            },
            'add_taxo_parameters': {
                'tax_lineages_file': os.path.join(temp_dir, "nonexistent_tax")
            },
        }
        
        with pytest.raises(FileNotFoundError, match="Required database files not found"):
            validate_database_paths(config)
