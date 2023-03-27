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

Essentially you can create a `ContinentContouring` object using rotation files, the continent features (polygons), a contour resolution, a gap threshold, an area threshold and an age range. And then ask it to reconstruct the polygons to an `age` and contour them into continents:

```
continent_contouring = ContinentContouring(...)
...
contoured_continents = continent_contouring.get_contoured_continents(age)
```

Where the returned `contoured_continents` is a list of `ContouredContinent`. These are all the contoured continents on the globe at the specified age. Each of these can be used to query a continent perimeter, area and whether arbitrary points are contained inside it. If you want to query distance to the continent you can first retrieve its polygons and then query distance to those (each contoured continent is actually one or more pyGPlates polygons representing its boundary between land and ocean, eg, an exterior polygon with interior holes).

*A note on the input parameters mentioned above*... The contour resolution determines how finely tessellated the contour outlines are. The gap threshold controls how close the reconstructed polygons can be to each other to be joined together – this helps remove narrow channels/gaps, however it also has the effect of expanding the continent outwards. The area threshold excludes any polygon boundary (of a contoured continent) with area below the threshold - remember that a contoured continent can contain more than one boundary - so for example if a hole inside a continent has an area below the threshold then the hole disappears.

> Documentation for other modules is coming soon...
