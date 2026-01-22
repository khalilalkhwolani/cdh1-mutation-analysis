"""
CDH1 Mutation Analysis Package

A comprehensive bioinformatics pipeline for analyzing CDH1 (E-cadherin) 
protein mutations across species using sequence alignment, phylogenetic 
analysis, and deep learning prediction.
"""

__version__ = "1.0.0"
__author__ = "CDH1 Research Team"
__email__ = "research@cdh1analysis.org"

from .pipeline import CDH1Pipeline
from .utils.config import Config
from .utils.logger import setup_logger

__all__ = [
    "CDH1Pipeline",
    "Config", 
    "setup_logger"
]