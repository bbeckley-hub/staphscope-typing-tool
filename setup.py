from setuptools import setup, find_packages
import os
from setuptools.command.install import install

# Load long description from README.md
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)
        # Set execute permissions on MLST binary
        import stat
        mlst_path = os.path.join(self.install_lib, 'staphscope', 'tools', 'mlst', 'bin', 'mlst')
        if os.path.exists(mlst_path):
            os.chmod(mlst_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

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
            "tools/mlst/**/*",
            "tools/spatyper/**/*",
            "tools/sccmecfinder/**/*",
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
    cmdclass={
        'install': PostInstallCommand,
    },
)
