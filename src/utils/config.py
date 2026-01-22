"""Configuration management for CDH1 analysis pipeline."""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class Config:
    """Configuration manager for the CDH1 analysis pipeline."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration.
        
        Args:
            config_path: Path to YAML configuration file
        """
        self.config_path = config_path or "config/default.yaml"
        self._config = self._load_config()
        self._validate_config()
        
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as file:
                config = yaml.safe_load(file)
                logger.info(f"Configuration loaded from {self.config_path}")
                return config
        except FileNotFoundError:
            logger.error(f"Configuration file not found: {self.config_path}")
            raise
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML configuration: {e}")
            raise
            
    def _validate_config(self) -> None:
        """Validate configuration structure and required fields."""
        required_sections = ['paths', 'species', 'alignment', 'deep_learning']
        
        for section in required_sections:
            if section not in self._config:
                raise ValueError(f"Missing required configuration section: {section}")
                
        # Validate paths exist or can be created
        for path_key, path_value in self._config['paths'].items():
            path = Path(path_value)
            if not path.exists():
                path.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created directory: {path}")
                
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value using dot notation.
        
        Args:
            key: Configuration key (e.g., 'paths.data_root')
            default: Default value if key not found
            
        Returns:
            Configuration value
        """
        keys = key.split('.')
        value = self._config
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            if default is not None:
                return default
            raise KeyError(f"Configuration key not found: {key}")
            
    def set(self, key: str, value: Any) -> None:
        """
        Set configuration value using dot notation.
        
        Args:
            key: Configuration key (e.g., 'paths.data_root')
            value: Value to set
        """
        keys = key.split('.')
        config = self._config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
            
        config[keys[-1]] = value
        
    def save(self, output_path: Optional[str] = None) -> None:
        """
        Save configuration to YAML file.
        
        Args:
            output_path: Output file path (defaults to original path)
        """
        output_path = output_path or self.config_path
        
        with open(output_path, 'w', encoding='utf-8') as file:
            yaml.dump(self._config, file, default_flow_style=False, indent=2)
            
        logger.info(f"Configuration saved to {output_path}")
        
    @property
    def species_list(self) -> list:
        """Get list of species names."""
        return list(self._config['species'].keys())
        
    @property
    def species_files(self) -> Dict[str, str]:
        """Get mapping of species to file names."""
        return {
            species: info['file'] 
            for species, info in self._config['species'].items()
        }
        
    def get_species_info(self, species: str) -> Dict[str, Any]:
        """
        Get complete information for a species.
        
        Args:
            species: Species name
            
        Returns:
            Species information dictionary
        """
        if species not in self._config['species']:
            raise ValueError(f"Unknown species: {species}")
            
        return self._config['species'][species]
        
    def __getitem__(self, key: str) -> Any:
        """Allow dictionary-style access."""
        return self.get(key)
        
    def __setitem__(self, key: str, value: Any) -> None:
        """Allow dictionary-style assignment."""
        self.set(key, value)
        
    def __contains__(self, key: str) -> bool:
        """Check if key exists in configuration."""
        try:
            self.get(key)
            return True
        except KeyError:
            return False