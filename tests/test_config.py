"""Tests for configuration management."""

import pytest
import tempfile
import yaml
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from src.utils.config import Config


class TestConfig:
    """Test configuration management functionality."""
    
    def test_config_initialization(self):
        """Test config initialization with default file."""
        # Create temporary config file
        config_data = {
            'paths': {'data_root': 'data/'},
            'species': {'human': {'name': 'Homo sapiens', 'file': 'human.fasta'}},
            'alignment': {'pairwise': {'method': 'globalxx'}},
            'deep_learning': {'model': {'type': 'LSTM'}}
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name
        
        try:
            config = Config(config_path)
            assert config.get('paths.data_root') == 'data/'
            assert config.get('species.human.name') == 'Homo sapiens'
        finally:
            Path(config_path).unlink()
    
    def test_config_get_set(self):
        """Test getting and setting configuration values."""
        config_data = {
            'test': {'nested': {'value': 42}},
            'simple': 'test_value'
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name
        
        try:
            config = Config(config_path)
            
            # Test getting values
            assert config.get('test.nested.value') == 42
            assert config.get('simple') == 'test_value'
            assert config.get('nonexistent', 'default') == 'default'
            
            # Test setting values
            config.set('test.nested.new_value', 100)
            assert config.get('test.nested.new_value') == 100
            
            # Test dictionary-style access
            assert config['simple'] == 'test_value'
            config['new_key'] = 'new_value'
            assert config['new_key'] == 'new_value'
            
        finally:
            Path(config_path).unlink()
    
    def test_species_methods(self):
        """Test species-specific methods."""
        config_data = {
            'species': {
                'human': {'name': 'Homo sapiens', 'file': 'human.fasta'},
                'mouse': {'name': 'Mus musculus', 'file': 'mouse.fasta'}
            }
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name
        
        try:
            config = Config(config_path)
            
            # Test species list
            assert set(config.species_list) == {'human', 'mouse'}
            
            # Test species files
            files = config.species_files
            assert files['human'] == 'human.fasta'
            assert files['mouse'] == 'mouse.fasta'
            
            # Test species info
            human_info = config.get_species_info('human')
            assert human_info['name'] == 'Homo sapiens'
            assert human_info['file'] == 'human.fasta'
            
        finally:
            Path(config_path).unlink()
    
    def test_config_validation(self):
        """Test configuration validation."""
        # Missing required sections
        invalid_config = {'paths': {}}
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(invalid_config, f)
            config_path = f.name
        
        try:
            with pytest.raises(ValueError, match="Missing required configuration section"):
                Config(config_path)
        finally:
            Path(config_path).unlink()
    
    def test_config_file_not_found(self):
        """Test handling of missing config file."""
        with pytest.raises(FileNotFoundError):
            Config("nonexistent_config.yaml")
    
    def test_invalid_yaml(self):
        """Test handling of invalid YAML."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("invalid: yaml: content: [")
            config_path = f.name
        
        try:
            with pytest.raises(yaml.YAMLError):
                Config(config_path)
        finally:
            Path(config_path).unlink()
    
    def test_key_error_handling(self):
        """Test handling of missing keys."""
        config_data = {'existing': 'value'}
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name
        
        try:
            config = Config(config_path)
            
            # Should raise KeyError for missing key without default
            with pytest.raises(KeyError):
                config.get('nonexistent.key')
            
            # Should return default for missing key with default
            assert config.get('nonexistent.key', 'default') == 'default'
            
        finally:
            Path(config_path).unlink()
    
    def test_contains_method(self):
        """Test __contains__ method."""
        config_data = {'existing': {'nested': 'value'}}
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = f.name
        
        try:
            config = Config(config_path)
            
            assert 'existing' in config
            assert 'existing.nested' in config
            assert 'nonexistent' not in config
            assert 'existing.nonexistent' not in config
            
        finally:
            Path(config_path).unlink()