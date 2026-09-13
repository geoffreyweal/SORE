# Installation: Setting Up SORE and Pre-Requisites Packages

In this article, we will look at how to install the SORE program and all the requisites required for this program.

## Pre-requisites

### Python 3 and pip3

This program is written in Python 3. You will need Python 3.4 or later. See [the Python website](https://www.python.org/downloads/) if you do not have Python 3 installed.

You will also need ``pip3`` to install SORE and its requisites. See [Installing pip](https://pip.pypa.io/en/stable/installation/) if you do not have it.

### Requisite Python packages

SORE requires the following python packages. Installing SORE with ``pip3`` (below) will install all of them for you, so you do not normally need to install them by hand.

* ``numpy``
* ``ase`` (version 3.19.0 or later)
* ``packaging``
* ``tqdm``
* ``xlsxwriter``

### SUMELF

The ``SUMELF`` program contains the supporting methods that are shared across the programs in the grand scheme, and **every program in the grand scheme depends on it**. SORE is no exception.

``SUMELF`` is not published on PyPI, so install it from GitHub:

```bash
pip3 install --upgrade --user git+https://github.com/geoffreyweal/SUMELF.git
```

See the [SUMELF Installation webpage](https://geoffreyweal.github.io/SUMELF/Installation) for more information.

!!! note

	SORE depends only on ``SUMELF``. It does not require any of the other programs in the grand scheme to be installed. Installing SORE with the command below pulls ``SUMELF`` in automatically.

## Setting up the SORE Program

SORE is not published on PyPI or conda, so install it from GitHub:

```bash
pip3 install --upgrade --user git+https://github.com/geoffreyweal/SORE.git
```

To check that SORE installed correctly, type the following into a python session:

```python
from SORE import Run_SORE
```

This should import without an error. Note that SORE is a python library rather than a terminal command, so there is no ``sore`` command to run.


## Upgrading SORE

To upgrade to the latest version of SORE, run the install command again:

```bash
pip3 install --upgrade --user git+https://github.com/geoffreyweal/SORE.git
```
