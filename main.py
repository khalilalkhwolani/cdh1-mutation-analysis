#!/usr/bin/env python3
"""
CDH1 Mutation Analysis Pipeline - Main Entry Point

This script provides the main entry point for running the CDH1 mutation analysis pipeline.
It supports various analysis modes and configuration options.

Usage:
    python main.py --config config/default.yaml
    python main.py --mode alignment --species human,chimp
    python main.py --mode full --output results/
"""

import argparse
import sys
from pathlib import Path
import traceback

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.pipeline import CDH1Pipeline
from src.utils.config import Config
from src.utils.logger import setup_logger


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="CDH1 Mutation Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline with default configuration
  python main.py
  
  # Use custom configuration
  python main.py --config config/production.yaml
  
  # Run only alignment analysis
  python main.py --mode alignment
  
  # Analyze specific species
  python main.py --species human,chimp,mouse
  
  # Set custom output directory
  python main.py --output results/experiment_1/
  
  # Enable debug logging
  python main.py --log-level DEBUG
  
  # Generate report only (requires existing results)
  python main.py --mode report
        """
    )
    
    parser.add_argument(
        "--config", "-c",
        type=str,
        default="config/default.yaml",
        help="Path to configuration YAML file (default: config/default.yaml)"
    )
    
    parser.add_argument(
        "--mode", "-m",
        type=str,
        choices=["full", "alignment", "phylogenetic", "deeplearning", "report"],
        default="full",
        help="Analysis mode to run (default: full)"
    )
    
    parser.add_argument(
        "--species", "-s",
        type=str,
        help="Comma-separated list of species to analyze (e.g., human,chimp,mouse)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output directory for results (overrides config)"
    )
    
    parser.add_argument(
        "--log-level", "-l",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level (default: INFO)"
    )
    
    parser.add_argument(
        "--log-file",
        type=str,
        help="Log file path (auto-generated if not specified)"
    )
    
    parser.add_argument(
        "--no-console-log",
        action="store_true",
        help="Disable console logging"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing"
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing results"
    )
    
    parser.add_argument(
        "--version", "-v",
        action="version",
        version="CDH1 Mutation Analysis Pipeline v1.0.0"
    )
    
    return parser.parse_args()


def validate_arguments(args):
    """Validate command line arguments."""
    errors = []
    
    # Check if config file exists
    if not Path(args.config).exists():
        errors.append(f"Configuration file not found: {args.config}")
    
    # Validate species list if provided
    if args.species:
        valid_species = ["human", "chimp", "mouse", "rat"]
        species_list = [s.strip() for s in args.species.split(",")]
        invalid_species = [s for s in species_list if s not in valid_species]
        if invalid_species:
            errors.append(f"Invalid species: {invalid_species}. Valid options: {valid_species}")
    
    # Check output directory permissions if specified
    if args.output:
        output_path = Path(args.output)
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            errors.append(f"Cannot create output directory: {args.output}")
    
    return errors


def main():
    """Main function."""
    # Parse arguments
    args = parse_arguments()
    
    # Validate arguments
    validation_errors = validate_arguments(args)
    if validation_errors:
        print("Error: Invalid arguments:")
        for error in validation_errors:
            print(f"  - {error}")
        sys.exit(1)
    
    # Setup logging
    logger = setup_logger(
        level=args.log_level,
        log_file=args.log_file,
        console_output=not args.no_console_log
    )
    
    logger.info("=" * 60)
    logger.info("CDH1 Mutation Analysis Pipeline Starting")
    logger.info("=" * 60)
    logger.info(f"Mode: {args.mode}")
    logger.info(f"Configuration: {args.config}")
    
    try:
        # Load configuration
        config = Config(args.config)
        
        # Override configuration with command line arguments
        if args.output:
            config.set('paths.results', args.output)
            logger.info(f"Output directory set to: {args.output}")
        
        if args.species:
            species_list = [s.strip() for s in args.species.split(",")]
            # Filter species configuration to only include specified species
            all_species = config.get('species')
            filtered_species = {k: v for k, v in all_species.items() if k in species_list}
            config.set('species', filtered_species)
            logger.info(f"Analyzing species: {species_list}")
        
        # Initialize pipeline
        pipeline = CDH1Pipeline(config, dry_run=args.dry_run, force=args.force)
        
        # Run analysis based on mode
        if args.mode == "full":
            logger.info("Running full analysis pipeline")
            results = pipeline.run_full_analysis()
            
        elif args.mode == "alignment":
            logger.info("Running alignment analysis only")
            results = pipeline.run_alignment_analysis()
            
        elif args.mode == "phylogenetic":
            logger.info("Running phylogenetic analysis only")
            results = pipeline.run_phylogenetic_analysis()
            
        elif args.mode == "deeplearning":
            logger.info("Running deep learning analysis only")
            results = pipeline.run_deeplearning_analysis()
            
        elif args.mode == "report":
            logger.info("Generating report from existing results")
            results = pipeline.generate_report()
            
        else:
            raise ValueError(f"Unknown mode: {args.mode}")
        
        # Print summary
        if results and not args.dry_run:
            logger.info("=" * 60)
            logger.info("Analysis completed successfully!")
            logger.info("=" * 60)
            
            if 'summary' in results:
                summary = results['summary']
                logger.info("Results Summary:")
                for key, value in summary.items():
                    logger.info(f"  {key}: {value}")
            
            if 'output_files' in results:
                logger.info("Output files generated:")
                for file_path in results['output_files']:
                    logger.info(f"  - {file_path}")
        
        elif args.dry_run:
            logger.info("Dry run completed - no files were modified")
        
    except KeyboardInterrupt:
        logger.warning("Analysis interrupted by user")
        sys.exit(1)
        
    except Exception as e:
        logger.error(f"Analysis failed with error: {e}")
        logger.debug("Full traceback:")
        logger.debug(traceback.format_exc())
        sys.exit(1)
    
    logger.info("Pipeline execution completed")


if __name__ == "__main__":
    main()