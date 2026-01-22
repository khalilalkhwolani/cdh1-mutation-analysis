"""Setup script for CDH1 Mutation Analysis package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

# Read requirements
requirements_path = Path(__file__).parent / "requirements.txt"
requirements = []
if requirements_path.exists():
    with open(requirements_path, 'r', encoding='utf-8') as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name="cdh1-mutation-analysis",
    version="1.0.0",
    author="CDH1 Research Team",
    author_email="research@cdh1analysis.org",
    description="Bioinformatics pipeline for CDH1 protein mutation analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-username/CDH1-Mutation-Analysis",
    project_urls={
        "Bug Reports": "https://github.com/your-username/CDH1-Mutation-Analysis/issues",
        "Source": "https://github.com/your-username/CDH1-Mutation-Analysis",
        "Documentation": "https://cdh1-mutation-analysis.readthedocs.io/",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "flake8>=6.0.0",
            "black>=23.7.0",
            "mypy>=1.5.1",
            "pre-commit>=3.3.3",
        ],
        "docs": [
            "sphinx>=7.1.2",
            "sphinx-rtd-theme>=1.3.0",
            "mkdocs>=1.5.2",
            "mkdocs-material>=9.1.21",
        ],
        "jupyter": [
            "jupyter>=1.0.0",
            "jupyterlab>=4.0.5",
            "ipykernel>=6.25.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "cdh1-analysis=main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yaml", "*.yml", "*.json", "*.md"],
    },
    zip_safe=False,
    keywords=[
        "bioinformatics",
        "cdh1",
        "e-cadherin", 
        "mutation-analysis",
        "sequence-alignment",
        "phylogenetics",
        "deep-learning",
        "protein-analysis",
        "gastric-cancer",
    ],
)