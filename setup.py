import os
from setuptools import setup, find_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "taxanalyser",
            "taxanalyser.VERSION",
        )
    ) as f:
        return f.readline().strip()
    

def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="taxanalyser",
    packages=find_packages(),
    url="https://github.com/gbouras13/taxanalyser",
    python_requires=">=3.9",
    description="Snakemake and Snaketool pipeline to taxonomically profile ONT long read metagenomics with sourmash",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="George Bouras",
    author_email="george.bouras@adelaide.edu.au",
    data_files=get_data_files(),
    py_modules=["taxanalyser"],
    install_requires=[
        "snaketool-utils>=0.0.2",
        "snakemake>=7.14.0",
        "pyyaml>=6.0",
        "Click>=8.1.3",
        "metasnek>=0.0.6",
        "attrmap>=0.0.5"
    ],
    entry_points={
        "console_scripts": [
            "taxanalyser=taxanalyser.__main__:main"
        ]
    },
    include_package_data=True,
)
