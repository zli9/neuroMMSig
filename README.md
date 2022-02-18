<!--

<p align="center">
  <img src="https://github.com///raw/main/docs/source/logo.png" height="150">
</p>

-->

<h1 align="center">
  mechanrich
</h1>


<p align="center">
    <a href="https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig/actions?query=workflow%3ATests">
        <img alt="Tests" src="https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig/workflows/Tests/badge.svg" />
    </a>
    <a href="https://test.pypi.org/project/mechanrich">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/enrichment" />
    </a>
    <a href="https://test.pypi.org/project/mechanrich">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/mechanrich" />
    </a>
    <a href="https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/mechanrich" />
    </a>
    <a href='https://enrichment.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/mechanrich/badge/?version=latest' alt='Documentation Status' />
    </a>
    <a href="https://codecov.io/gh///branch/main">
        <img src="https://codecov.io/gh///branch/main/graph/badge.svg" alt="Codecov status" />
    </a>  
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img alt="Cookiecutter template from @cthoyt" src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" /> 
    </a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' />
    </a>
    <a href="https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig/blob/main/.github/CODE_OF_CONDUCT.md">
        <img src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" alt="Contributor Covenant"/>
    </a>
</p>

Python package for mechanism enrichment using reverse causal reasoning (RCR).

---
## 💪 Getting Started

### Prerequisites

- Gene expression experiment data
- Pathway data of interest
- Mapping table containing regulation information

### 1. Get enriched network

```python
from mechanrich import mechanrich as mec

graph = mec.Graph()

# get an enriched pathway network
graph.plot_full_network(output="pathway_network.pdf", dpi=300)

# get a HYP network of SMAD2 gene
graph.plot_hyp_network(gene="SMAD2", output="hyp_network_SMAD2.pdf", dpi=300)
```

### 2. Get statistics of RCR inference

```python
stat = mec.RCRstat()

# get a statistics table
stat.get_stat(output="stat.txt")

# get concordance of HYP network of SMAD2 gene
conc = stat.get_gene_conc(gene="SMAD2")

# get richness of HYP network of SMAD2 gene
rich = stat.get_gene_rich(gene="SMAD2")
```

### Command Line Interface

The mechanrich command line tool is automatically installed. It can
be used from the shell with the `--help` flag to show all subcommands:

```shell
$ mechanrich --help
```

## 🚀 Installation

The most recent release can be installed from
[TestPyPI](https://test.pypi.org/project/mechanrich/) with:

```bash
$ pip install --index-url https://test.pypi.org/simple/ --extra-index-url  https://pypi.org/simple/ mechanrich==0.0.3.dev0
```

The most recent code and data can be installed directly from GitHub with:

```bash
$ pip install https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig.git
```

## 📦 Project Structure

```angular2html
│  .bumpversion.cfg
│  .gitignore
│  .readthedocs.yml
│  LICENSE
│  MANIFEST.in
│  pyproject.toml
│  README.md
│  setup.cfg
│  tox.ini
│  
├─.github
│          
├─data
│      
├─docs
│          
├─src
│  │  __init__.py
│  │  
│  └─mechanrich
│          __init__.py
│          __main__.py
│          cli.py
│          constants.py
│          mechanrich.py
│          preprocessing.py
│          reader.py
│          startup.py
│          utils.py
│          
└─tests
    │  __init__.py
    │  test_mechanrich.py
    │  
    └─test_data
```

## 👋 Attribution

### ⚖️ License

The code in this package is licensed under the MIT License.


### 📖 Citation

Catlett, N.L., Bargnesi, A.J., Ungerer, S. et al. Reverse causal reasoning: applying qualitative causal knowledge to the interpretation of high-throughput data. BMC Bioinformatics 14, 340 (2013). https://doi.org/10.1186/1471-2105-14-340

<!--

### 🎁 Support

This project has been supported by the following organizations (in alphabetical order):

- [Harvard Program in Therapeutic Science - Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/)

-->

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>


The final section of the README is for if you want to get involved by making a code contribution.

### Development Installation

To install in development mode, use the following:

```bash
$ git clone https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig.git
$ cd 
$ pip install -e .
```

### 🥼 Testing

After cloning the repository and installing `tox` with `pip install tox`, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
$ tox
```

Additionally, these tests are automatically re-run with each commit in a [GitHub Action](https://github.com///actions?query=workflow%3ATests).

### 📖 Building the Documentation

The documentation can be built locally using the following:

```shell
$ git clone https://github.com/zli9/Mechanism-enrichment-using-NeuroMMSig.git
$ cd 
$ tox -e docs
$ open docs/build/html/index.html
```

The documentation automatically installs the package as well as the `docs`
extra specified in the [`setup.cfg`](setup.cfg). `sphinx` plugins
like `texext` can be added there. Additionally, they need to be added to the
`extensions` list in [`docs/source/conf.py`](docs/source/conf.py).

### 📦 Making a Release

After installing the package in development mode and installing
`tox` with `pip install tox`, the commands for making a new release are contained within the `finish` environment
in `tox.ini`. Run the following from the shell:

```shell
$ tox -e finish
```

This script does the following:

1. Uses [Bump2Version](https://github.com/c4urself/bump2version) to switch the version number in the `setup.cfg`,
   `src/mechanrich/version.py`, and [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using [`build`](https://github.com/pypa/build)
3. Uploads to PyPI using [`twine`](https://github.com/pypa/twine). Be sure to have a `.pypirc` file configured to avoid the need for manual input at this
   step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use `tox -e bumpversion minor` after.
   </details>
