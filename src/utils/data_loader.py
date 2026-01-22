"""Data loading utilities for CDH1 analysis pipeline."""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import LoggerMixin
from .config import Config


class DataLoader(LoggerMixin):
    """Data loader for sequence files and analysis results."""
    
    def __init__(self, config: Config):
        """
        Initialize data loader.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.sequences_path = Path(config.get('paths.sequences'))
        
    def load_sequence(self, species: str) -> SeqRecord:
        """
        Load sequence for a specific species.
        
        Args:
            species: Species name
            
        Returns:
            BioPython SeqRecord object
            
        Raises:
            FileNotFoundError: If sequence file not found
            ValueError: If species not configured
        """
        if species not in self.config.species_list:
            raise ValueError(f"Unknown species: {species}")
            
        species_info = self.config.get_species_info(species)
        file_path = self.sequences_path / species_info['file']
        
        if not file_path.exists():
            raise FileNotFoundError(f"Sequence file not found: {file_path}")
            
        try:
            sequence = SeqIO.read(file_path, "fasta")
            self.logger.info(f"Loaded sequence for {species}: {len(sequence)} amino acids")
            return sequence
        except Exception as e:
            self.logger.error(f"Error loading sequence for {species}: {e}")
            raise
            
    def load_all_sequences(self) -> Dict[str, SeqRecord]:
        """
        Load sequences for all configured species.
        
        Returns:
            Dictionary mapping species names to SeqRecord objects
        """
        sequences = {}
        
        for species in self.config.species_list:
            try:
                sequences[species] = self.load_sequence(species)
            except Exception as e:
                self.logger.warning(f"Failed to load sequence for {species}: {e}")
                
        self.logger.info(f"Loaded sequences for {len(sequences)} species")
        return sequences
        
    def validate_sequences(self, sequences: Dict[str, SeqRecord]) -> bool:
        """
        Validate loaded sequences.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            True if all sequences are valid
        """
        valid = True
        
        for species, seq in sequences.items():
            # Check sequence length
            if len(seq) == 0:
                self.logger.error(f"Empty sequence for {species}")
                valid = False
                continue
                
            # Check for valid amino acids
            valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
            seq_aa = set(str(seq.seq).upper())
            invalid_aa = seq_aa - valid_aa
            
            if invalid_aa:
                self.logger.warning(
                    f"Invalid amino acids in {species}: {invalid_aa}"
                )
                
            self.logger.info(
                f"{species}: {len(seq)} amino acids, "
                f"composition: {dict(sorted([(aa, str(seq.seq).count(aa)) for aa in valid_aa if str(seq.seq).count(aa) > 0]))}"
            )
                
        return valid
        
    def save_sequences(self, sequences: Dict[str, SeqRecord], output_dir: Optional[Path] = None) -> None:
        """
        Save sequences to FASTA files.
        
        Args:
            sequences: Dictionary of sequences
            output_dir: Output directory (defaults to sequences path)
        """
        output_dir = output_dir or self.sequences_path
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for species, seq in sequences.items():
            species_info = self.config.get_species_info(species)
            output_file = output_dir / species_info['file']
            
            SeqIO.write(seq, output_file, "fasta")
            self.logger.info(f"Saved sequence for {species} to {output_file}")
            
    def load_alignment_results(self, results_file: str) -> pd.DataFrame:
        """
        Load alignment results from CSV file.
        
        Args:
            results_file: Path to results CSV file
            
        Returns:
            DataFrame with alignment results
        """
        results_path = Path(self.config.get('paths.results')) / results_file
        
        if not results_path.exists():
            raise FileNotFoundError(f"Results file not found: {results_path}")
            
        try:
            df = pd.read_csv(results_path)
            self.logger.info(f"Loaded alignment results from {results_path}")
            return df
        except Exception as e:
            self.logger.error(f"Error loading alignment results: {e}")
            raise
            
    def save_results(self, data: pd.DataFrame, filename: str, output_dir: Optional[Path] = None) -> None:
        """
        Save results to CSV file.
        
        Args:
            data: DataFrame to save
            filename: Output filename
            output_dir: Output directory (defaults to results path)
        """
        output_dir = output_dir or Path(self.config.get('paths.results'))
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_dir / filename
        data.to_csv(output_file, index=False)
        self.logger.info(f"Saved results to {output_file}")
        
    def get_file_info(self) -> Dict[str, Dict[str, any]]:
        """
        Get information about available data files.
        
        Returns:
            Dictionary with file information
        """
        info = {
            'sequences': {},
            'results': {},
            'models': {}
        }
        
        # Sequence files
        for species in self.config.species_list:
            species_info = self.config.get_species_info(species)
            file_path = self.sequences_path / species_info['file']
            
            if file_path.exists():
                info['sequences'][species] = {
                    'file': species_info['file'],
                    'path': str(file_path),
                    'size': file_path.stat().st_size,
                    'exists': True
                }
            else:
                info['sequences'][species] = {
                    'file': species_info['file'],
                    'path': str(file_path),
                    'exists': False
                }
                
        # Results files
        results_dir = Path(self.config.get('paths.results'))
        if results_dir.exists():
            for file_path in results_dir.glob('*.csv'):
                info['results'][file_path.stem] = {
                    'file': file_path.name,
                    'path': str(file_path),
                    'size': file_path.stat().st_size
                }
                
        # Model files
        models_dir = Path(self.config.get('paths.models'))
        if models_dir.exists():
            for file_path in models_dir.glob('*.keras'):
                info['models'][file_path.stem] = {
                    'file': file_path.name,
                    'path': str(file_path),
                    'size': file_path.stat().st_size
                }
                
        return info