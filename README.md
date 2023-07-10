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

The following Python packages are required:

- [`numpy`](http://numpy.org)
- [`scipy`](https://scipy.org)
- [`pygplates`](http://gplates.org/docs/pygplates/pygplates_getting_started.html#installation)

__Optional dependencies__ for running the Notebooks:

- [`matplotlib`](https://matplotlib.org/)
- [`cartopy`](https://scitools.org.uk/cartopy/docs/latest/)

### Installing using pip

You can install `ptt` using the
[`pip package manager`](https://pypi.org/project/pip/) with either version of Python:

```bash
python2 -m pip install PlateTectonicTools
python3 -m pip install PlateTectonicTools
```

To install the latest version from GitHub, use:

```bash
pip3 install --no-cache-dir --upgrade git+https://github.com/EarthByte/PlateTectonicTools
pip install --no-cache-dir --upgrade git+https://github.com/EarthByte/PlateTectonicTools
```

### Installing using Docker

> Coming soon...


### API Documentation

#### continent_contours.py

The current contouring algorithm implements a landmass flood-fill to find all contours for a particular continuous landmass followed by the Marching Squares algorithm to create the contours for that landmass. This is done on the 3D globe and solves the following problems:

*	A single very large landmass can have more than one contour (polygon).
*	In some cases (particularly for 0-130Ma) the polygons are so large that the inside/outside region of each contour polygon becomes swapped.
*	Clipping at the dateline (with a 2D approach) causes perimeters to be larger than they should be.

Essentially you can create a `ContinentContouring` object using rotation files, some continent/craton features (polygons), a contour resolution to determine how finely tessellated the contour outlines should be, a buffer/gap threshold to expand continents outward, and two area thresholds to separately exclude small continental islands and small oceanic islands. It will then reconstruct the polygons to an `age` and contour them into continents. This will output the continent contours (as *polyline* continental-oceanic boundaries) and/or continent masks (a 2D NumPy boolean array of continental crust at each age):

```
from ptt.continent_contours import ContinentContouring

continent_contouring = ContinentContouring(...)
...
contoured_continents = continent_contouring.get_contoured_continents(age)
continent_mask = continent_contouring.get_continent_mask(age)
```

The returned `contoured_continents` is a list of `ContouredContinent`. And `continent_mask` is a 2D NumPy boolean array of continental crust (that can be converted to floating-point `0.0` and `1.0` values using `continent_mask.astype('float')`). The contoured continents can be used to query a continent perimeter, a continent area and whether arbitrary points are contained inside a continent. If you want to query distance to the continent-ocean boundaries you can first retrieve the contours (which are *polylines*) and then query distance to those.

> Documentation for other modules is coming soon...
