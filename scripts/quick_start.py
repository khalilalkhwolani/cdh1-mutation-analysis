#!/usr/bin/env python3
"""
Quick start script for CDH1 Mutation Analysis.

This script demonstrates basic usage of the pipeline with sample data.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from src.utils.config import Config
from src.utils.logger import setup_logger
from src.pipeline import CDH1Pipeline


def create_sample_data():
    """Create sample FASTA files for testing."""
    sequences = {
        'human_CDH1.fasta': '>sp|P12830|CADH1_HUMAN\nMCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY',
        'chimp_CDH1.fasta': '>tr|A0A2J8Q404|A0A2J8Q404_PANTR\nMCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWH',
        'mouse_CDH1.fasta': '>sp|P09803|CADH1_MOUSE\nMCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY',
        'rat_CDH1.fasta': '>sp|Q9R0T4|CADH1_RAT\nMCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY'
    }
    
    data_dir = Path('data/sequences')
    data_dir.mkdir(parents=True, exist_ok=True)
    
    for filename, content in sequences.items():
        with open(data_dir / filename, 'w') as f:
            f.write(content)
    
    print(f"Sample data created in {data_dir}")


def main():
    """Run quick start demonstration."""
    print("CDH1 Mutation Analysis - Quick Start")
    print("=" * 40)
    
    # Setup logging
    logger = setup_logger(level="INFO")
    
    # Create sample data
    create_sample_data()
    
    # Load configuration
    config = Config('config/default.yaml')
    
    # Initialize pipeline
    pipeline = CDH1Pipeline(config)
    
    # Run alignment analysis
    logger.info("Running alignment analysis...")
    results = pipeline.run_alignment_analysis()
    
    # Print results summary
    print("\nResults Summary:")
    print(pipeline.get_results_summary())
    
    print("\nQuick start completed successfully!")


if __name__ == "__main__":
    main()