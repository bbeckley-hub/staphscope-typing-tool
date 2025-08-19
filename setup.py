from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="staphscope-typing-tool",
    version="0.1.0",
    description="Staphscope Typing Tool — MLST + spa + SCCmec typing for Staphylococcus aureus",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/staphscope-typing-tool",
    author="Beckley Brown",
    author_email="brownbeckley94@gmail.com",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="bioinformatics, microbiology, staphylococcus, typing, mlst, spa, sccmec",
    packages=find_packages(),
    python_requires=">=3.6, <4",
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
    ],
    entry_points={
        "console_scripts": [
            "staphscope=staphscope.cli:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/bbeckley-hub/staphscope-typing-tool/issues",
        "Source": "https://github.com/bbeckley-hub/staphscope-typing-tool",
    },
)
