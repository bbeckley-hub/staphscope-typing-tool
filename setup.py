from setuptools import setup, find_packages
import os

# Read the README for long description
try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "Advanced Staphylococcus aureus Typing & Lineage Analysis Platform"

setup(
    name="staphscope",
    version="1.0.0",
    author="Brown Beckley",
    author_email="brownbeckley94@gmail.com",
    description="Advanced Staphylococcus aureus Typing & Lineage Analysis Platform",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/staphscope-typing-tool",
    packages=find_packages(),
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
        "Programming Language :: Python :: 3.12",
    ],
    python_requires="==3.10.*",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.80",  
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
        "click>=8.0.0"
    ],
    entry_points={
        "console_scripts": [
            "staphscope=staphscope.staphscope:main",
        ],
    },
    include_package_data=True,
    package_data={
        'staphscope': [
            '**/*',
        ]
    },
    zip_safe=False,
)
