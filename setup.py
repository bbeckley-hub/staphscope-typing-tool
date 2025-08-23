from setuptools import setup, find_packages

setup(
    name="staphscope",
    version="0.2.0",
    description="A typing tool for Staphylococcus aureus (spa, MLST, SCCmec typing)",
    author="Beckley Brown",
    author_email="brownbeckley94@gmail.com",
    packages=find_packages(),
    package_data={
        'staphscope': [
            'database/mlst/**/*',
            'database/spatyper/**/*',
            'database/sccmecfinder/**/*',
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
    ],
    python_requires=">=3.7",
)
