# Contributing to StaphScope

Thank you for your interest in contributing to StaphScope! We welcome contributions from researchers, clinicians, and developers working to combat antimicrobial resistance. Whether you're fixing bugs, adding features, improving documentation, or reporting issues, your help is invaluable.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Features](#suggesting-features)
- [Contributing Code](#contributing-code)
- [Development Setup](#development-setup)
- [Coding Style](#coding-style)
- [Pull Request Process](#pull-request-process)
- [Questions?](#questions)

## Code of Conduct

This project adheres to a [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to brownbeckley94@gmail.com.

## How Can I Contribute?

### Reporting Bugs
If you find a bug, please create an issue with:
- **Clear title** describing the problem
- **Steps to reproduce** with example commands
- **Expected vs actual behavior**
- **Environment details**: OS, StaphScope version, Conda version
- **Error logs** or screenshots

Before submitting, please check if the issue already exists.

### Suggesting Features
We're always looking to improve StaphScope. Feature suggestions should include:
- **Clear description** of the feature
- **Use case** - why it would be valuable
- **Potential implementation** ideas (if you have any)
- **Alternative solutions** you've considered

### Contributing Code

#### Development Setup

1. **Fork and clone the repository:**
```bash
git clone https://github.com/YOUR_USERNAME/staphscope-typing-tool.git
cd staphscope-typing-tool
```

2. **Set up the development environment:**
```bash
# Using Conda (recommended)
conda env create -f environment.yml
conda activate staphscope
pip install -e .
```

3. **Install additional databases (for testing):**
```bash
# AMR database
staphscope --force-update-amr-db

# ABRicate databases
abricate --setupdb
```

#### Coding Style

- Follow **PEP 8** for Python code
- Use **meaningful variable names** (e.g., `sample_names` not `sn`)
- Add **docstrings** to functions and classes
- Keep functions **modular and focused** on single tasks
- Add **comments** for complex logic
- Use **type hints** where possible

Example:
```python
def process_sample(sample_path: str, output_dir: str, threads: int = 4) -> dict:
    """
    Process a single genome sample through all StaphScope modules.
    
    Args:
        sample_path: Path to the input FASTA file
        output_dir: Directory to write results
        threads: Number of CPU threads to use
    
    Returns:
        Dictionary containing all analysis results
    """
    # Implementation here
```

#### Commit Messages

Use clear, descriptive commit messages:
- **Good:** `Add dynamic grouping by MLST in Ultimate Reporter`
- **Good:** `Fix memory leak in AMR module`
- **Avoid:** `Update code` or `Fix bug`

#### Testing

Before submitting, ensure your changes:
- Don't break existing functionality
- Include tests if adding new features
- Pass all existing tests

## Pull Request Process

1. **Create a feature branch:**
```bash
git checkout -b feature/your-feature-name
```

2. **Make your changes** and commit them with clear messages.

3. **Push to your fork:**
```bash
git push origin feature/your-feature-name
```

4. **Open a Pull Request** against the `main` branch.

5. **PR Description should include:**
   - What the changes do
   - Why they're needed
   - Any breaking changes
   - Screenshots for UI changes (if applicable)

6. **Review Process:**
   - Maintainers will review your PR
   - Address any feedback or questions
   - PR will be merged once approved

## Questions?

If you have questions about contributing:
- Open a **GitHub Discussion**
- Email Brown Beckley at **brownbeckley94@gmail.com**
- Check the **README** for more information

## Recognition

All contributors will be acknowledged in the project. We appreciate your help in making StaphScope better!

---

**Thank you for contributing to the fight against antimicrobial resistance!** 🔬🦠

