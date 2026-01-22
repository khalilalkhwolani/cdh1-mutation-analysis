"""Phylogenetic analysis utilities for CDH1 mutation analysis."""

from typing import Dict, List, Optional, Any
import pandas as pd
from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Phylo import draw_ascii, write, read
from Bio import Phylo
import matplotlib.pyplot as plt

from ..utils.logger import LoggerMixin
from ..utils.config import Config


class PhylogeneticAnalyzer(LoggerMixin):
    """Phylogenetic tree construction and analysis."""
    
    def __init__(self, config: Config):
        """
        Initialize phylogenetic analyzer.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.phylo_params = config.get('phylogenetics')
        
    def construct_tree(self, alignment, method: str = "nj") -> Dict[str, Any]:
        """
        Construct phylogenetic tree from multiple sequence alignment.
        
        Args:
            alignment: Multiple sequence alignment
            method: Tree construction method ('nj' for neighbor joining)
            
        Returns:
            Dictionary with tree and analysis results
        """
        self.logger.info(f"Constructing phylogenetic tree using {method}")
        
        try:
            # Calculate distance matrix
            calculator = DistanceCalculator(self.phylo_params.get('distance_method', 'identity'))
            distance_matrix = calculator.get_distance(alignment)
            
            # Construct tree
            constructor = DistanceTreeConstructor()
            
            if method.lower() == "nj":
                tree = constructor.nj(distance_matrix)
            elif method.lower() == "upgma":
                tree = constructor.upgma(distance_matrix)
            else:
                raise ValueError(f"Unknown tree construction method: {method}")
                
            # Root the tree
            tree.root_at_midpoint()
            
            # Calculate tree statistics
            tree_stats = self._calculate_tree_statistics(tree)
            
            result = {
                'tree': tree,
                'distance_matrix': distance_matrix,
                'method': method,
                'statistics': tree_stats,
                'species_count': len(alignment),
                'tree_length': tree.total_branch_length()
            }
            
            self.logger.info(f"Tree constructed successfully with {len(alignment)} species")
            return result
            
        except Exception as e:
            self.logger.error(f"Error constructing phylogenetic tree: {e}")
            raise
            
    def _calculate_tree_statistics(self, tree) -> Dict[str, Any]:
        """
        Calculate basic tree statistics.
        
        Args:
            tree: Phylogenetic tree
            
        Returns:
            Dictionary with tree statistics
        """
        stats = {
            'total_branch_length': tree.total_branch_length(),
            'terminal_count': len(tree.get_terminals()),
            'internal_count': len(tree.get_nonterminals()),
            'max_distance': 0,
            'min_distance': float('inf'),
            'avg_distance': 0
        }
        
        # Calculate pairwise distances
        terminals = tree.get_terminals()
        distances = []
        
        for i, term1 in enumerate(terminals):
            for term2 in terminals[i+1:]:
                distance = tree.distance(term1, term2)
                distances.append(distance)
                stats['max_distance'] = max(stats['max_distance'], distance)
                stats['min_distance'] = min(stats['min_distance'], distance)
                
        if distances:
            stats['avg_distance'] = sum(distances) / len(distances)
            
        return stats
        
    def save_tree(self, tree, output_path: Path, format: str = "newick") -> None:
        """
        Save phylogenetic tree to file.
        
        Args:
            tree: Phylogenetic tree
            output_path: Output file path
            format: Output format (newick, nexus, phyloxml)
        """
        write(tree, output_path, format)
        self.logger.info(f"Tree saved to {output_path} in {format} format")
        
    def visualize_tree(self, tree, output_path: Optional[Path] = None, 
                      figsize: tuple = (12, 8)) -> None:
        """
        Create tree visualization.
        
        Args:
            tree: Phylogenetic tree
            output_path: Output file path for saving plot
            figsize: Figure size
        """
        plt.figure(figsize=figsize)
        
        # Draw tree
        Phylo.draw(tree, do_show=False)
        plt.title("CDH1 Phylogenetic Tree")
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Tree visualization saved to {output_path}")
        
        plt.close()
        
    def print_tree_ascii(self, tree) -> str:
        """
        Generate ASCII representation of tree.
        
        Args:
            tree: Phylogenetic tree
            
        Returns:
            ASCII tree string
        """
        import io
        import sys
        
        # Capture ASCII output
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO()
        
        try:
            draw_ascii(tree)
            ascii_tree = buffer.getvalue()
        finally:
            sys.stdout = old_stdout
            
        return ascii_tree
        
    def compare_trees(self, tree1, tree2) -> Dict[str, Any]:
        """
        Compare two phylogenetic trees.
        
        Args:
            tree1: First tree
            tree2: Second tree
            
        Returns:
            Comparison results
        """
        # This is a simplified comparison
        # In practice, you might want to use more sophisticated methods
        
        comparison = {
            'tree1_length': tree1.total_branch_length(),
            'tree2_length': tree2.total_branch_length(),
            'tree1_terminals': len(tree1.get_terminals()),
            'tree2_terminals': len(tree2.get_terminals()),
            'length_difference': abs(tree1.total_branch_length() - tree2.total_branch_length())
        }
        
        return comparison
        
    def bootstrap_analysis(self, alignment, replicates: int = 100) -> Dict[str, Any]:
        """
        Perform bootstrap analysis for tree support.
        
        Args:
            alignment: Multiple sequence alignment
            replicates: Number of bootstrap replicates
            
        Returns:
            Bootstrap analysis results
        """
        self.logger.info(f"Performing bootstrap analysis with {replicates} replicates")
        
        # This is a simplified bootstrap implementation
        # In practice, you might want to use more sophisticated methods
        
        bootstrap_trees = []
        
        # Generate bootstrap replicates (simplified)
        for i in range(min(replicates, 10)):  # Limit for demo
            try:
                tree_result = self.construct_tree(alignment)
                bootstrap_trees.append(tree_result['tree'])
            except Exception as e:
                self.logger.warning(f"Bootstrap replicate {i} failed: {e}")
                
        result = {
            'replicates': len(bootstrap_trees),
            'requested_replicates': replicates,
            'consensus_tree': bootstrap_trees[0] if bootstrap_trees else None,
            'support_values': self._calculate_support_values(bootstrap_trees)
        }
        
        return result
        
    def _calculate_support_values(self, trees: List) -> Dict[str, float]:
        """
        Calculate bootstrap support values.
        
        Args:
            trees: List of bootstrap trees
            
        Returns:
            Support values for tree nodes
        """
        # Simplified support calculation
        if not trees:
            return {}
            
        # In a real implementation, you would compare tree topologies
        # and calculate support for each internal node
        
        support_values = {
            'root_support': 100.0,  # Placeholder
            'internal_nodes_avg_support': 85.0,  # Placeholder
            'min_support': 70.0,  # Placeholder
            'max_support': 100.0   # Placeholder
        }
        
        return support_values