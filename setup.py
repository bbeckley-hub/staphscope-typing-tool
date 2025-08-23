from setuptools import setup, find_packages

# Load long description from README.md
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="staphscope",
    version="0.2.1",
    description="Unified MLST + spa + SCCmec typing tool for Staphylococcus aureus",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Brown Beckley",
    author_email="brownbeckley94@gmail.com",
    url="https://github.com/bbeckley-hub/staphscope-typing-tool",
    packages=find_packages(),
    package_data={
        "staphscope": [
            "database/mlst/**/*",
            "database/spatyper/**/*",
            "database/sccmecfinder/**/*",
        ]
    },
    include_package_data=True,
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "click>=8.0.0",
    ],
    entry_points={
        "console_scripts": [
            "staphscope=staphscope.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    keywords="staphylococcus aureus typing mlst spa sccmec bioinformatics",
    project_urls={
        "Source": "https://github.com/bbeckley-hub/staphscope-typing-tool",
        "Bug Reports": "https://github.com/bbeckley-hub/staphscope-typing-tool/issues",
    },
)
