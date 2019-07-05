# PlateTectonicTools

Python tools for plate tectonic research.

This repository contains the Python package `ptt` (short for Plate Tectonic Tools) which provides a collection of common plate tectonic functionality that researchers can use in their workflows. It is primarily built on top of the [pyGPlates](http://gplates.org/docs/pygplates/index.html) Python library.

There are also some Jupyter notebooks to demonstrate usage of the **ptt** package:

* Calculating average seafloor spreading rates of mid-ocean ridges
* Calculating subducting rate of paleo rasters along subduction zones
* Predicting CO2 from age and bottom water temperature.

*Note:* This repository is only in its initial stages. More functionality and examples will be provided in future.


## Documentation

There are some demonstrations of `ptt` that may be installed from the package itself by running:

```python
import ptt
ptt.install_documentation(path="PTT-Notebooks")
```

## Installation

### Dependencies

You will need **Python 2.7**.
Also, the following packages are required:

- [`numpy`](http://numpy.org)
- [`scipy`](https://scipy.org)
- [`pygplates`](http://gplates.org/docs/pygplates/pygplates_getting_started.html#installation)

__Optional dependencies__ for running the Notebooks:

- [`matplotlib`](https://matplotlib.org/)

### Installing using pip

You can install `ptt` using the
[`pip package manager`](https://pypi.org/project/pip/) with either version of Python:

```bash
python2 -m pip install PlateTectonicTools
python3 -m pip install PlateTectonicTools
```

### Installing using Docker

> Coming soon...


### API Documentation

> Coming soon...
