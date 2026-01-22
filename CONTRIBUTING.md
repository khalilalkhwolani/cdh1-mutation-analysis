# Contributing to CDH1 Mutation Analysis

We welcome contributions to the CDH1 Mutation Analysis project! This document provides guidelines for contributing to the project.

## 🚀 Getting Started

### Prerequisites
- Python 3.8 or higher
- Git
- Basic knowledge of bioinformatics and Python

### Development Setup

1. **Fork the repository**
   ```bash
   git clone https://github.com/your-username/CDH1-Mutation-Analysis.git
   cd CDH1-Mutation-Analysis
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install development dependencies**
   ```bash
   pip install -r requirements-dev.txt
   ```

4. **Install pre-commit hooks**
   ```bash
   pre-commit install
   ```

5. **Run tests to ensure everything works**
   ```bash
   pytest tests/
   ```

## 📝 How to Contribute

### Reporting Bugs

Before creating bug reports, please check the existing issues to avoid duplicates. When creating a bug report, include:

- **Clear description** of the problem
- **Steps to reproduce** the issue
- **Expected vs actual behavior**
- **Environment details** (OS, Python version, etc.)
- **Error messages** and stack traces
- **Sample data** if applicable

### Suggesting Enhancements

Enhancement suggestions are welcome! Please include:

- **Clear description** of the enhancement
- **Use case** and motivation
- **Proposed implementation** (if you have ideas)
- **Potential impact** on existing functionality

### Pull Requests

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Follow the coding standards (see below)
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes**
   ```bash
   # Run tests
   pytest tests/
   
   # Run linting
   flake8 src/
   black src/
   
   # Type checking
   mypy src/
   ```

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "feat: add new alignment algorithm"
   ```

5. **Push and create pull request**
   ```bash
   git push origin feature/your-feature-name
   ```

## 🎯 Coding Standards

### Python Style Guide

We follow PEP 8 with some modifications:

- **Line length**: 88 characters (Black default)
- **Imports**: Use absolute imports, group by standard library, third-party, local
- **Docstrings**: Use Google style docstrings
- **Type hints**: Required for all public functions

### Code Formatting

We use automated code formatting:

```bash
# Format code
black src/ tests/

# Sort imports
isort src/ tests/

# Lint code
flake8 src/ tests/

# Type checking
mypy src/
```

### Documentation

- **Docstrings**: All public functions, classes, and modules must have docstrings
- **Comments**: Use comments to explain complex logic, not obvious code
- **README updates**: Update README.md for significant changes
- **API docs**: Update API documentation for new features

### Testing

- **Unit tests**: Write tests for all new functions
- **Integration tests**: Add tests for new pipeline components
- **Test coverage**: Aim for >90% test coverage
- **Test data**: Use small, synthetic test datasets

Example test structure:
```python
def test_pairwise_alignment():
    """Test pairwise alignment functionality."""
    # Arrange
    seq1 = create_test_sequence("ACDEFG")
    seq2 = create_test_sequence("ACDEFH")
    aligner = PairwiseAligner(test_config)
    
    # Act
    result = aligner.align_sequences(seq1, seq2, "test1", "test2")
    
    # Assert
    assert result['similarity_percent'] > 80
    assert result['mutations_count'] == 1
```

## 🧬 Bioinformatics Guidelines

### Sequence Data
- **File formats**: Support standard formats (FASTA, FASTQ, etc.)
- **Validation**: Always validate input sequences
- **Error handling**: Handle malformed sequences gracefully
- **Memory efficiency**: Consider memory usage for large sequences

### Algorithms
- **Performance**: Profile algorithms with realistic data sizes
- **Accuracy**: Validate against known benchmarks
- **Parameters**: Make algorithm parameters configurable
- **Documentation**: Document algorithm choices and limitations

### Results
- **Reproducibility**: Ensure results are reproducible with same inputs
- **Validation**: Include statistical validation where appropriate
- **Formats**: Use standard output formats (CSV, JSON, etc.)
- **Metadata**: Include analysis metadata in results

## 📊 Project Structure

When adding new features, follow the project structure:

```
src/
├── bioinformatics/     # Sequence analysis algorithms
├── deeplearning/       # ML/DL models and training
├── utils/             # Common utilities
└── visualization/     # Plotting and visualization

tests/
├── unit/              # Unit tests
├── integration/       # Integration tests
└── fixtures/          # Test data and fixtures

docs/
├── api/               # API documentation
├── tutorials/         # Usage tutorials
└── methodology.md     # Scientific methodology
```

## 🔄 Development Workflow

### Branch Naming
- `feature/description` - New features
- `bugfix/description` - Bug fixes
- `docs/description` - Documentation updates
- `refactor/description` - Code refactoring

### Commit Messages
Follow conventional commits format:

```
type(scope): description

[optional body]

[optional footer]
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

Examples:
```
feat(alignment): add MUSCLE alignment algorithm
fix(pipeline): handle empty sequence files
docs(readme): update installation instructions
test(alignment): add pairwise alignment tests
```

### Code Review Process

1. **Self-review**: Review your own code before submitting
2. **Automated checks**: Ensure all CI checks pass
3. **Peer review**: At least one reviewer approval required
4. **Testing**: Verify tests pass and coverage is maintained
5. **Documentation**: Ensure documentation is updated

## 🧪 Testing Guidelines

### Test Categories

1. **Unit Tests**: Test individual functions and classes
2. **Integration Tests**: Test component interactions
3. **End-to-end Tests**: Test complete workflows
4. **Performance Tests**: Test with realistic data sizes

### Test Data

- Use small, synthetic datasets for unit tests
- Include edge cases (empty sequences, single amino acid, etc.)
- Test with real data samples (anonymized)
- Document test data sources and characteristics

### Running Tests

```bash
# All tests
pytest

# Specific test file
pytest tests/test_alignment.py

# With coverage
pytest --cov=src --cov-report=html

# Performance tests
pytest tests/performance/ -m slow
```

## 📚 Documentation

### Types of Documentation

1. **Code documentation**: Docstrings and comments
2. **API documentation**: Auto-generated from docstrings
3. **User guides**: How to use the software
4. **Developer guides**: How to contribute and extend
5. **Scientific documentation**: Methodology and validation

### Documentation Standards

- **Clarity**: Write for your target audience
- **Completeness**: Cover all public APIs
- **Examples**: Include code examples
- **Updates**: Keep documentation in sync with code

## 🤝 Community Guidelines

### Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help newcomers learn and contribute
- Acknowledge contributions from others

### Communication

- **Issues**: Use GitHub issues for bugs and feature requests
- **Discussions**: Use GitHub discussions for questions and ideas
- **Email**: Contact maintainers for sensitive issues

## 🏆 Recognition

Contributors will be recognized in:

- **README.md**: Contributors section
- **CHANGELOG.md**: Release notes
- **Documentation**: Author credits
- **Publications**: Co-authorship for significant contributions

## 📞 Getting Help

If you need help:

1. **Check documentation**: README, docs/, and code comments
2. **Search issues**: Look for similar problems
3. **Ask questions**: Create a GitHub discussion
4. **Contact maintainers**: Email for urgent issues

## 🎉 Thank You!

Thank you for contributing to the CDH1 Mutation Analysis project! Your contributions help advance bioinformatics research and improve tools for the scientific community.

---

For questions about contributing, please create a GitHub discussion or contact the maintainers.