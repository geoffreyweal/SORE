# Issues and Troubleshooting

This page lists problems that can come up when running the SORE program, and what to do about them.

## ``ModuleNotFoundError: No module named 'SUMELF'``

SORE depends on the ``SUMELF`` program. Install it into the same python environment:

```bash
pip3 install --upgrade --user git+https://github.com/geoffreyweal/SUMELF.git
```

See [Installation](Installation.md) for more information.

If SUMELF *is* installed and you still see this, check that you are not running from the folder that holds your repository clones — see [the section below](#importerror-or-modulenotfounderror-when-working-in-the-folder-that-holds-your-clones).

## ``Error: SORE_settings['mode'] must be set to either: ...``

``SORE_settings`` must contain a ``mode`` entry, set to one of ``'standard'``, ``'analytical_energetic_disorder'`` or ``'numeric_energetic_disorder'``. See [How To Use The SORE Program](Using_The_SORE_Program.md) for what each mode does.

```python
SORE_settings = {'mode': 'standard'}
```

## There is no ``sore`` command

There is not meant to be. SORE is a python library, not a terminal program: you run it by calling ``Run_SORE`` from a ``Run_SORE.py`` script. See [How To Use The SORE Program](Using_The_SORE_Program.md).

## Setting up a crystal is very slow

Building the neighbourhood of every molecule in the crystal is the slow part of the calculation. Raise ``no_of_cpus_for_setup``:

```python
Run_SORE(EKMC_settings, SORE_settings=SORE_settings, no_of_cpus_for_setup=8)
```

If you want to try several ``SORE_settings`` on the same crystal, set ``save_initial_data=True`` on the first run so that the neighbourhood and coupling data can be reused rather than recomputed each time.

## SORE and EKMC give different diffusion coefficients

This is expected to some degree: SORE evaluates a sum-over-rates equation while EKMC averages explicit simulated trajectories. Large disagreement is worth investigating, and usually means the two runs were not given the same ``EKMC_settings``, or that the kinetic model or coupling settings differ between them.


## ``ImportError`` or ``ModuleNotFoundError`` when working in the folder that holds your clones

If you keep your clones of these programs together in one folder, like this:

~~~
Solar_Cell_Project/
|-- ACSD/
|-- ECCP/
|-- SUMELF/
|-- ...
~~~

then running a command or a python script **from that folder** fails with an error such as:

~~~
ImportError: cannot import name '__version__' from 'SORE' (unknown location)
~~~

or:

~~~
ModuleNotFoundError: No module named 'SUMELF'
~~~

even though you installed everything correctly.

This happens because python searches the current directory first. The ``SORE`` folder sitting there is the **git repository**, not the python package: the package is one level further in, at ``SORE/SORE``. The repository folder has no ``__init__.py``, so python treats it as an empty namespace package, finds no code inside it, and stops looking. Your installed copy is never reached.

The fix is to work from anywhere other than that folder. Change into the directory holding the crystals you are working on and run from there, which is what you would normally be doing anyway:

~~~bash
cd /path/to/my_crystals
python3 Run_SORE.py
~~~

!!! tip

	Only the folder that *directly* contains the repository folders is affected. Subfolders of it are fine, and so is any unrelated directory.

## Something else has gone wrong

If you have found a problem that is not covered here, please [open an issue on GitHub](https://github.com/geoffreyweal/SORE/issues) describing:

* what you were trying to do,
* the settings you used,
* and the full error message you received.
