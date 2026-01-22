"""Sequence alignment utilities for CDH1 analysis."""

from typing import Dict, List, Tuple, Optional
import pandas as pd
from pathlib import Path
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess
import tempfile
import os

from ..utils.logger import LoggerMixin
from ..utils.config import Config


class PairwiseAligner(LoggerMixin):
    """Pairwise sequence alignment using BioPython."""
    
    def __init__(self, config: Config):
        """
        Initialize pairwise aligner.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.alignment_params = config.get('alignment.pairwise')
        
    def align_sequences(
        self, 
        seq1: SeqRecord, 
        seq2: SeqRecord,
        species1: str,
        species2: str
    ) -> Dict[str, any]:
        """
        Perform pairwise alignment between two sequences.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            species1: Name of first species
            species2: Name of second species
            
        Returns:
            Dictionary with alignment results
        """
        self.logger.info(f"Aligning {species1} vs {species2}")
        
        try:
            # Perform global alignment
            alignments = pairwise2.align.globalxx(
                seq1.seq, 
                seq2.seq,
                match=self.alignment_params.get('match_score', 2),
                mismatch=self.alignment_params.get('mismatch_score', -1),
                gap_open=self.alignment_params.get('gap_open', -2),
                gap_extend=self.alignment_params.get('gap_extend', -0.5)
            )
            
            if not alignments:
                raise ValueError("No alignments found")
                
            best_alignment = alignments[0]
            
            # Calculate similarity
            matches = sum(
                a == b for a, b in zip(best_alignment.seqA, best_alignment.seqB)
                if a != '-' and b != '-'
            )
            
            # Calculate alignment length (excluding gaps)
            alignment_length = sum(
                1 for a, b in zip(best_alignment.seqA, best_alignment.seqB)
                if a != '-' and b != '-'
            )
            
            similarity = (matches / alignment_length) * 100 if alignment_length > 0 else 0
            
            # Detect mutations
            mutations = self._detect_mutations(best_alignment)
            
            result = {
                'species1': species1,
                'species2': species2,
                'similarity_percent': round(similarity, 2),
                'matches': matches,
                'alignment_length': alignment_length,
                'mutations_count': len(mutations),
                'mutations': mutations[:10],  # First 10 mutations
                'alignment_score': best_alignment.score,
                'seq1_length': len(seq1.seq),
                'seq2_length': len(seq2.seq),
                'aligned_seq1': str(best_alignment.seqA),
                'aligned_seq2': str(best_alignment.seqB)
            }
            
            self.logger.info(
                f"{species1} vs {species2}: {similarity:.2f}% similarity, "
                f"{len(mutations)} mutations"
            )
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error in pairwise alignment: {e}")
            raise
            
    def _detect_mutations(self, alignment) -> List[Tuple[int, str, str]]:
        """
        Detect mutations from alignment.
        
        Args:
            alignment: BioPython alignment object
            
        Returns:
            List of mutations as (position, amino_acid1, amino_acid2)
        """
        mutations = []
        
        for i, (a, b) in enumerate(zip(alignment.seqA, alignment.seqB)):
            if a != b and a != "-" and b != "-":
                mutations.append((i + 1, a, b))
                
        return mutations
        
    def align_all_pairs(self, sequences: Dict[str, SeqRecord]) -> pd.DataFrame:
        """
        Perform pairwise alignment for all species pairs.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            DataFrame with all pairwise alignment results
        """
        results = []
        species_list = list(sequences.keys())
        
        for i, species1 in enumerate(species_list):
            for j, species2 in enumerate(species_list):
                if i < j:  # Avoid duplicate comparisons
                    try:
                        result = self.align_sequences(
                            sequences[species1],
                            sequences[species2],
                            species1,
                            species2
                        )
                        results.append({
                            'comparison': f"{species1} vs {species2}",
                            'species1': species1,
                            'species2': species2,
                            'similarity_percent': result['similarity_percent'],
                            'mutations_count': result['mutations_count'],
                            'alignment_score': result['alignment_score']
                        })
                    except Exception as e:
                        self.logger.error(
                            f"Failed alignment {species1} vs {species2}: {e}"
                        )
                        
        return pd.DataFrame(results)


class MultipleAligner(LoggerMixin):
    """Multiple sequence alignment using CLUSTAL Omega."""
    
    def __init__(self, config: Config):
        """
        Initialize multiple aligner.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.alignment_params = config.get('alignment.multiple')
        
    def align_sequences(self, sequences: Dict[str, SeqRecord]) -> Dict[str, any]:
        """
        Perform multiple sequence alignment.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            Dictionary with alignment results
        """
        self.logger.info(f"Performing multiple alignment for {len(sequences)} sequences")
        
        try:
            # Create temporary input file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
                for species, seq in sequences.items():
                    # Use species info for proper sequence ID
                    species_info = self.config.get_species_info(species)
                    seq_id = species_info.get('uniprot_id', species)
                    temp_input.write(f">{seq_id}\n{seq.seq}\n")
                temp_input_path = temp_input.name
                
            # Create temporary output file
            temp_output = tempfile.NamedTemporaryFile(suffix='.aln', delete=False)
            temp_output_path = temp_output.name
            temp_output.close()
            
            # Run CLUSTAL Omega (if available) or use BioPython fallback
            try:
                self._run_clustal_omega(temp_input_path, temp_output_path)
                alignment = AlignIO.read(temp_output_path, "clustal")
            except (FileNotFoundError, subprocess.CalledProcessError):
                self.logger.warning("CLUSTAL Omega not available, using BioPython fallback")
                alignment = self._fallback_alignment(sequences)
                
            # Calculate distance matrix
            distance_matrix = self._calculate_distance_matrix(alignment)
            identity_matrix = self._calculate_identity_matrix(alignment)
            
            result = {
                'alignment': alignment,
                'distance_matrix': distance_matrix,
                'identity_matrix': identity_matrix,
                'sequence_count': len(alignment),
                'alignment_length': alignment.get_alignment_length(),
                'species_names': [record.id for record in alignment]
            }
            
            # Clean up temporary files
            os.unlink(temp_input_path)
            if os.path.exists(temp_output_path):
                os.unlink(temp_output_path)
                
            self.logger.info(
                f"Multiple alignment completed: {len(alignment)} sequences, "
                f"{alignment.get_alignment_length()} positions"
            )
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error in multiple alignment: {e}")
            raise
            
    def _run_clustal_omega(self, input_file: str, output_file: str) -> None:
        """
        Run CLUSTAL Omega command line tool.
        
        Args:
            input_file: Input FASTA file
            output_file: Output alignment file
        """
        clustalo_cline = ClustalOmegaCommandline(
            infile=input_file,
            outfile=output_file,
            verbose=True,
            auto=True,
            outfmt="clustal"
        )
        
        stdout, stderr = clustalo_cline()
        
        if stderr:
            self.logger.warning(f"CLUSTAL Omega stderr: {stderr}")
            
    def _fallback_alignment(self, sequences: Dict[str, SeqRecord]):
        """
        Fallback alignment method using BioPython.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            Multiple sequence alignment
        """
        # Simple progressive alignment fallback
        # This is a basic implementation - in production, consider using
        # more sophisticated methods or ensuring CLUSTAL Omega is available
        
        from Bio.Align import MultipleSeqAlignment
        from Bio.SeqRecord import SeqRecord
        
        # Convert sequences to list
        seq_list = []
        for species, seq in sequences.items():
            species_info = self.config.get_species_info(species)
            seq_id = species_info.get('uniprot_id', species)
            new_record = SeqRecord(seq.seq, id=seq_id, description=species_info['name'])
            seq_list.append(new_record)
            
        # Create alignment (this is a placeholder - real implementation would
        # need proper multiple alignment algorithm)
        alignment = MultipleSeqAlignment(seq_list)
        
        self.logger.warning("Using fallback alignment - results may not be optimal")
        return alignment
        
    def _calculate_distance_matrix(self, alignment) -> pd.DataFrame:
        """
        Calculate distance matrix from alignment.
        
        Args:
            alignment: Multiple sequence alignment
            
        Returns:
            Distance matrix as DataFrame
        """
        from Bio.Phylo.TreeConstruction import DistanceCalculator
        
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        
        # Convert to DataFrame
        names = distance_matrix.names
        data = []
        
        for i, name1 in enumerate(names):
            row = []
            for j, name2 in enumerate(names):
                if j <= i:
                    row.append(distance_matrix[name1, name2])
                else:
                    row.append(distance_matrix[name2, name1])
            data.append(row)
            
        df = pd.DataFrame(data, index=names, columns=names)
        return df
        
    def _calculate_identity_matrix(self, alignment) -> pd.DataFrame:
        """
        Calculate identity matrix from alignment.
        
        Args:
            alignment: Multiple sequence alignment
            
        Returns:
            Identity matrix as DataFrame
        """
        names = [record.id for record in alignment]
        identity_data = []
        
        for i in range(len(alignment)):
            row = []
            for j in range(len(alignment)):
                identity = self._calculate_identity(alignment[i].seq, alignment[j].seq)
                row.append(round(identity, 2))
            identity_data.append(row)
            
        df = pd.DataFrame(identity_data, index=names, columns=names)
        return df
        
    def _calculate_identity(self, seq1, seq2) -> float:
        """
        Calculate percentage identity between two sequences.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Percentage identity
        """
        matches = 0
        length = 0
        
        for a, b in zip(seq1, seq2):
            if a != '-' and b != '-':  # Ignore gaps
                length += 1
                if a == b:
                    matches += 1
                    
        return (matches / length) * 100 if length > 0 else 0
        
    def save_alignment(self, alignment, output_path: Path, format: str = "clustal") -> None:
        """
        Save alignment to file.
        
        Args:
            alignment: Multiple sequence alignment
            output_path: Output file path
            format: Output format (clustal, fasta, phylip, etc.)
        """
        AlignIO.write(alignment, output_path, format)
        self.logger.info(f"Alignment saved to {output_path} in {format} format")