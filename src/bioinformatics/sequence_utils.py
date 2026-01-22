"""Sequence analysis utilities for CDH1 mutation analysis."""

from typing import Dict, List, Tuple, Optional
import pandas as pd
from collections import Counter
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from ..utils.logger import LoggerMixin
from ..utils.config import Config


class SequenceAnalyzer(LoggerMixin):
    """Comprehensive sequence analysis utilities."""
    
    def __init__(self, config: Config):
        """
        Initialize sequence analyzer.
        
        Args:
            config: Configuration object
        """
        self.config = config
        
    def analyze_sequence(self, sequence: SeqRecord, species: str) -> Dict[str, any]:
        """
        Perform comprehensive analysis of a single sequence.
        
        Args:
            sequence: Sequence to analyze
            species: Species name
            
        Returns:
            Dictionary with analysis results
        """
        self.logger.info(f"Analyzing sequence for {species}")
        
        seq_str = str(sequence.seq)
        protein_analysis = ProteinAnalysis(seq_str)
        
        try:
            result = {
                'species': species,
                'sequence_id': sequence.id,
                'description': sequence.description,
                'length': len(sequence),
                'molecular_weight': round(molecular_weight(sequence.seq, seq_type='protein'), 2),
                'amino_acid_composition': dict(protein_analysis.get_amino_acids_percent()),
                'secondary_structure': protein_analysis.secondary_structure_fraction(),
                'instability_index': round(protein_analysis.instability_index(), 2),
                'isoelectric_point': round(protein_analysis.isoelectric_point(), 2),
                'gravy': round(protein_analysis.gravy(), 3),  # Grand average of hydropathy
                'aromaticity': round(protein_analysis.aromaticity(), 3),
                'flexibility': protein_analysis.flexibility(),
                'charge_at_pH': {
                    'pH_7': round(protein_analysis.charge_at_pH(7.0), 2),
                    'pH_7.4': round(protein_analysis.charge_at_pH(7.4), 2)
                }
            }
            
            # Add domain analysis if available
            result['domains'] = self._predict_domains(seq_str)
            
            # Add conservation analysis
            result['conservation_score'] = self._calculate_conservation_score(seq_str)
            
            self.logger.info(
                f"{species}: {len(sequence)} aa, MW: {result['molecular_weight']} Da, "
                f"pI: {result['isoelectric_point']}"
            )
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error analyzing sequence for {species}: {e}")
            raise
            
    def analyze_all_sequences(self, sequences: Dict[str, SeqRecord]) -> pd.DataFrame:
        """
        Analyze all sequences and return summary DataFrame.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            DataFrame with analysis results
        """
        results = []
        
        for species, sequence in sequences.items():
            try:
                analysis = self.analyze_sequence(sequence, species)
                
                # Flatten the result for DataFrame
                flat_result = {
                    'species': species,
                    'length': analysis['length'],
                    'molecular_weight': analysis['molecular_weight'],
                    'isoelectric_point': analysis['isoelectric_point'],
                    'instability_index': analysis['instability_index'],
                    'gravy': analysis['gravy'],
                    'aromaticity': analysis['aromaticity'],
                    'charge_pH7': analysis['charge_at_pH']['pH_7'],
                    'charge_pH7.4': analysis['charge_at_pH']['pH_7.4'],
                    'helix_fraction': analysis['secondary_structure'][0],
                    'turn_fraction': analysis['secondary_structure'][1],
                    'sheet_fraction': analysis['secondary_structure'][2],
                    'conservation_score': analysis['conservation_score']
                }
                
                results.append(flat_result)
                
            except Exception as e:
                self.logger.warning(f"Failed to analyze {species}: {e}")
                
        return pd.DataFrame(results)
        
    def _predict_domains(self, sequence: str) -> List[Dict[str, any]]:
        """
        Predict protein domains (simplified implementation).
        
        Args:
            sequence: Protein sequence
            
        Returns:
            List of predicted domains
        """
        # This is a simplified domain prediction
        # In production, you might want to use tools like InterProScan, Pfam, etc.
        
        domains = []
        
        # E-cadherin specific domains (based on known structure)
        if len(sequence) > 700:  # Typical E-cadherin length
            domains.append({
                'name': 'Signal peptide',
                'start': 1,
                'end': 23,
                'confidence': 0.8
            })
            
            # Extracellular cadherin domains (EC1-EC5)
            ec_domains = [
                (24, 134, 'EC1'),
                (135, 245, 'EC2'), 
                (246, 356, 'EC3'),
                (357, 467, 'EC4'),
                (468, 578, 'EC5')
            ]
            
            for start, end, name in ec_domains:
                if end <= len(sequence):
                    domains.append({
                        'name': name,
                        'start': start,
                        'end': end,
                        'confidence': 0.9
                    })
                    
            # Transmembrane domain
            if len(sequence) > 600:
                domains.append({
                    'name': 'Transmembrane',
                    'start': 579,
                    'end': 601,
                    'confidence': 0.85
                })
                
            # Cytoplasmic domain
            if len(sequence) > 650:
                domains.append({
                    'name': 'Cytoplasmic',
                    'start': 602,
                    'end': len(sequence),
                    'confidence': 0.8
                })
                
        return domains
        
    def _calculate_conservation_score(self, sequence: str) -> float:
        """
        Calculate a simple conservation score based on amino acid properties.
        
        Args:
            sequence: Protein sequence
            
        Returns:
            Conservation score (0-1)
        """
        # This is a simplified conservation score
        # In practice, you'd compare against multiple homologs
        
        # Count conserved residues (based on typical E-cadherin conservation)
        conserved_positions = 0
        total_positions = len(sequence)
        
        # E-cadherin has highly conserved calcium-binding sites
        # and cell adhesion motifs - this is a simplified approximation
        
        # Look for conserved motifs
        conserved_motifs = [
            'DXD',  # Calcium binding
            'DXNDN',  # Calcium binding
            'HAV',  # Cell adhesion
            'LDRE'  # Cytoplasmic domain
        ]
        
        for motif in conserved_motifs:
            if motif in sequence:
                conserved_positions += len(motif)
                
        conservation_score = min(conserved_positions / total_positions, 1.0)
        return round(conservation_score, 3)
        
    def compare_sequences(self, sequences: Dict[str, SeqRecord]) -> Dict[str, any]:
        """
        Compare sequences and identify key differences.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            Comparison results
        """
        if len(sequences) < 2:
            raise ValueError("Need at least 2 sequences for comparison")
            
        # Get human sequence as reference (if available)
        reference_species = 'human' if 'human' in sequences else list(sequences.keys())[0]
        reference_seq = str(sequences[reference_species].seq)
        
        comparisons = {}
        
        for species, seq_record in sequences.items():
            if species == reference_species:
                continue
                
            seq_str = str(seq_record.seq)
            
            # Calculate differences
            differences = []
            min_length = min(len(reference_seq), len(seq_str))
            
            for i in range(min_length):
                if reference_seq[i] != seq_str[i]:
                    differences.append({
                        'position': i + 1,
                        'reference': reference_seq[i],
                        'variant': seq_str[i],
                        'type': self._classify_mutation(reference_seq[i], seq_str[i])
                    })
                    
            # Length differences
            length_diff = len(seq_str) - len(reference_seq)
            
            comparisons[species] = {
                'differences_count': len(differences),
                'differences': differences[:20],  # First 20 differences
                'length_difference': length_diff,
                'similarity_percent': ((min_length - len(differences)) / min_length) * 100,
                'reference_species': reference_species
            }
            
        return {
            'reference_species': reference_species,
            'comparisons': comparisons,
            'summary': self._generate_comparison_summary(comparisons)
        }
        
    def _classify_mutation(self, ref_aa: str, var_aa: str) -> str:
        """
        Classify mutation type based on amino acid properties.
        
        Args:
            ref_aa: Reference amino acid
            var_aa: Variant amino acid
            
        Returns:
            Mutation classification
        """
        # Amino acid property groups
        hydrophobic = set('AILMFPWV')
        polar = set('NQST')
        charged_positive = set('RHK')
        charged_negative = set('DE')
        aromatic = set('FWY')
        small = set('AGCS')
        
        ref_props = []
        var_props = []
        
        # Classify reference amino acid
        if ref_aa in hydrophobic:
            ref_props.append('hydrophobic')
        if ref_aa in polar:
            ref_props.append('polar')
        if ref_aa in charged_positive:
            ref_props.append('positive')
        if ref_aa in charged_negative:
            ref_props.append('negative')
        if ref_aa in aromatic:
            ref_props.append('aromatic')
        if ref_aa in small:
            ref_props.append('small')
            
        # Classify variant amino acid
        if var_aa in hydrophobic:
            var_props.append('hydrophobic')
        if var_aa in polar:
            var_props.append('polar')
        if var_aa in charged_positive:
            var_props.append('positive')
        if var_aa in charged_negative:
            var_props.append('negative')
        if var_aa in aromatic:
            var_props.append('aromatic')
        if var_aa in small:
            var_props.append('small')
            
        # Determine mutation type
        if set(ref_props) == set(var_props):
            return 'conservative'
        elif ('positive' in ref_props and 'negative' in var_props) or \
             ('negative' in ref_props and 'positive' in var_props):
            return 'charge_reversal'
        elif ('hydrophobic' in ref_props and 'polar' in var_props) or \
             ('polar' in ref_props and 'hydrophobic' in var_props):
            return 'polarity_change'
        else:
            return 'non_conservative'
            
    def _generate_comparison_summary(self, comparisons: Dict[str, any]) -> Dict[str, any]:
        """
        Generate summary statistics for sequence comparisons.
        
        Args:
            comparisons: Comparison results
            
        Returns:
            Summary statistics
        """
        total_differences = sum(comp['differences_count'] for comp in comparisons.values())
        avg_similarity = sum(comp['similarity_percent'] for comp in comparisons.values()) / len(comparisons)
        
        # Mutation type distribution
        mutation_types = []
        for comp in comparisons.values():
            for diff in comp['differences']:
                mutation_types.append(diff['type'])
                
        mutation_distribution = dict(Counter(mutation_types))
        
        return {
            'total_species_compared': len(comparisons),
            'total_differences': total_differences,
            'average_similarity': round(avg_similarity, 2),
            'mutation_type_distribution': mutation_distribution,
            'most_similar_species': max(comparisons.keys(), 
                                      key=lambda x: comparisons[x]['similarity_percent']),
            'least_similar_species': min(comparisons.keys(), 
                                       key=lambda x: comparisons[x]['similarity_percent'])
        }


class MutationDetector(LoggerMixin):
    """Specialized mutation detection and analysis."""
    
    def __init__(self, config: Config):
        """
        Initialize mutation detector.
        
        Args:
            config: Configuration object
        """
        self.config = config
        
    def detect_pathogenic_mutations(self, sequence: str, species: str) -> List[Dict[str, any]]:
        """
        Detect potentially pathogenic mutations based on known CDH1 variants.
        
        Args:
            sequence: Protein sequence
            species: Species name
            
        Returns:
            List of potentially pathogenic mutations
        """
        # Known pathogenic CDH1 mutations (simplified list)
        known_pathogenic = [
            {'position': 23, 'mutation': 'A23V', 'pathogenicity': 'high'},
            {'position': 63, 'mutation': 'E63K', 'pathogenicity': 'medium'},
            {'position': 134, 'mutation': 'C134Y', 'pathogenicity': 'high'},
            {'position': 283, 'mutation': 'R283Q', 'pathogenicity': 'medium'},
            {'position': 356, 'mutation': 'L356P', 'pathogenicity': 'high'}
        ]
        
        detected_mutations = []
        
        # This is a simplified implementation
        # In practice, you'd compare against comprehensive mutation databases
        
        for mutation_info in known_pathogenic:
            pos = mutation_info['position'] - 1  # Convert to 0-based indexing
            if pos < len(sequence):
                expected_aa = mutation_info['mutation'][0]
                observed_aa = sequence[pos]
                
                if observed_aa != expected_aa:
                    detected_mutations.append({
                        'position': mutation_info['position'],
                        'reference': expected_aa,
                        'observed': observed_aa,
                        'mutation_name': f"{expected_aa}{mutation_info['position']}{observed_aa}",
                        'pathogenicity': mutation_info['pathogenicity'],
                        'species': species
                    })
                    
        return detected_mutations
        
    def analyze_mutation_hotspots(self, sequences: Dict[str, SeqRecord]) -> Dict[str, any]:
        """
        Identify mutation hotspots across species.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            Hotspot analysis results
        """
        if len(sequences) < 2:
            return {}
            
        # Use human as reference if available
        reference_species = 'human' if 'human' in sequences else list(sequences.keys())[0]
        reference_seq = str(sequences[reference_species].seq)
        
        # Count mutations at each position
        position_mutations = {}
        
        for species, seq_record in sequences.items():
            if species == reference_species:
                continue
                
            seq_str = str(seq_record.seq)
            min_length = min(len(reference_seq), len(seq_str))
            
            for i in range(min_length):
                if reference_seq[i] != seq_str[i]:
                    pos = i + 1
                    if pos not in position_mutations:
                        position_mutations[pos] = []
                    position_mutations[pos].append({
                        'species': species,
                        'reference': reference_seq[i],
                        'variant': seq_str[i]
                    })
                    
        # Identify hotspots (positions with mutations in multiple species)
        hotspots = {
            pos: mutations for pos, mutations in position_mutations.items()
            if len(mutations) > 1
        }
        
        # Sort by number of species with mutations
        sorted_hotspots = sorted(
            hotspots.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )
        
        return {
            'total_variable_positions': len(position_mutations),
            'hotspot_positions': len(hotspots),
            'top_hotspots': sorted_hotspots[:10],
            'reference_species': reference_species,
            'mutation_frequency': {
                pos: len(mutations) for pos, mutations in position_mutations.items()
            }
        }