"""Tests for alignment functionality."""

import pytest
import tempfile
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from src.bioinformatics.alignment import PairwiseAligner
from src.utils.config import Config


class TestPairwiseAligner:
    """Test pairwise alignment functionality."""
    
    @pytest.fixture
    def test_config(self):
        """Create test configuration."""
        config_data = {
            'paths': {'data_root': 'data/', 'results': 'results/'},
            'species': {
                'test1': {'name': 'Test Species 1', 'file': 'test1.fasta'},
                'test2': {'name': 'Test Species 2', 'file': 'test2.fasta'}
            },
            'alignment': {
                'pairwise': {
                    'method': 'globalxx',
                    'match_score': 2,
                    'mismatch_score': -1,
                    'gap_open': -2,
                    'gap_extend': -0.5
                }
            },
            'deep_learning': {'model': {'type': 'LSTM'}}
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            import yaml
            yaml.dump(config_data, f)
            config_path = f.name
        
        config = Config(config_path)
        Path(config_path).unlink()
        return config
    
    @pytest.fixture
    def test_sequences(self):
        """Create test sequences."""
        seq1 = SeqRecord(
            Seq("ACDEFGHIKLMNPQRSTVWY"),
            id="test1",
            description="Test sequence 1"
        )
        seq2 = SeqRecord(
            Seq("ACDEFGHIKLMNPQRSTVWH"),  # Last amino acid different
            id="test2", 
            description="Test sequence 2"
        )
        return seq1, seq2
    
    def test_pairwise_alignment_identical(self, test_config):
        """Test alignment of identical sequences."""
        aligner = PairwiseAligner(test_config)
        
        seq1 = SeqRecord(Seq("ACDEFG"), id="seq1")
        seq2 = SeqRecord(Seq("ACDEFG"), id="seq2")
        
        result = aligner.align_sequences(seq1, seq2, "species1", "species2")
        
        assert result['similarity_percent'] == 100.0
        assert result['mutations_count'] == 0
        assert result['species1'] == "species1"
        assert result['species2'] == "species2"
    
    def test_pairwise_alignment_different(self, test_config, test_sequences):
        """Test alignment of different sequences."""
        aligner = PairwiseAligner(test_config)
        seq1, seq2 = test_sequences
        
        result = aligner.align_sequences(seq1, seq2, "test1", "test2")
        
        assert result['similarity_percent'] < 100.0
        assert result['mutations_count'] > 0
        assert result['alignment_length'] > 0
        assert 'mutations' in result
        assert len(result['mutations']) <= 10  # Should limit to first 10
    
    def test_detect_mutations(self, test_config):
        """Test mutation detection."""
        aligner = PairwiseAligner(test_config)
        
        # Create mock alignment object
        class MockAlignment:
            def __init__(self, seqA, seqB, score=100):
                self.seqA = seqA
                self.seqB = seqB
                self.score = score
        
        alignment = MockAlignment("ACDEFG", "ACDHFG")
        mutations = aligner._detect_mutations(alignment)
        
        assert len(mutations) == 1
        assert mutations[0] == (3, 'D', 'H')  # Position 3, D->H
    
    def test_align_all_pairs(self, test_config):
        """Test alignment of all sequence pairs."""
        aligner = PairwiseAligner(test_config)
        
        sequences = {
            'seq1': SeqRecord(Seq("ACDEFG"), id="seq1"),
            'seq2': SeqRecord(Seq("ACDEFH"), id="seq2"),
            'seq3': SeqRecord(Seq("ACDEFG"), id="seq3")
        }
        
        results_df = aligner.align_all_pairs(sequences)
        
        # Should have 3 comparisons: seq1-seq2, seq1-seq3, seq2-seq3
        assert len(results_df) == 3
        assert 'comparison' in results_df.columns
        assert 'similarity_percent' in results_df.columns
        assert 'mutations_count' in results_df.columns
    
    def test_empty_sequences(self, test_config):
        """Test handling of empty sequences."""
        aligner = PairwiseAligner(test_config)
        
        seq1 = SeqRecord(Seq(""), id="empty1")
        seq2 = SeqRecord(Seq("ACDEFG"), id="normal")
        
        # Should handle empty sequences gracefully
        result = aligner.align_sequences(seq1, seq2, "empty", "normal")
        assert 'similarity_percent' in result
        assert 'mutations_count' in result
    
    def test_single_amino_acid(self, test_config):
        """Test alignment of single amino acid sequences."""
        aligner = PairwiseAligner(test_config)
        
        seq1 = SeqRecord(Seq("A"), id="single1")
        seq2 = SeqRecord(Seq("A"), id="single2")
        
        result = aligner.align_sequences(seq1, seq2, "single1", "single2")
        
        assert result['similarity_percent'] == 100.0
        assert result['mutations_count'] == 0
    
    def test_very_different_sequences(self, test_config):
        """Test alignment of very different sequences."""
        aligner = PairwiseAligner(test_config)
        
        seq1 = SeqRecord(Seq("AAAAAAAAAA"), id="seq1")
        seq2 = SeqRecord(Seq("GGGGGGGGGG"), id="seq2")
        
        result = aligner.align_sequences(seq1, seq2, "seq1", "seq2")
        
        assert result['similarity_percent'] == 0.0
        assert result['mutations_count'] == 10
    
    def test_alignment_parameters(self, test_config):
        """Test that alignment parameters are used correctly."""
        aligner = PairwiseAligner(test_config)
        
        # Check that parameters are loaded from config
        params = aligner.alignment_params
        assert params['match_score'] == 2
        assert params['mismatch_score'] == -1
        assert params['gap_open'] == -2
        assert params['gap_extend'] == -0.5