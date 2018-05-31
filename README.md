[![Build Status](https://travis-ci.org/NuGrid/NuGridPy.svg?branch=master)](https://travis-ci.org/NuGrid/NuGridPy)
[![PyPI version](https://badge.fury.io/py/NuGridpy.svg)](https://badge.fury.io/py/NuGridpy)

# NuGridPy

[NuGridPy](https://nugrid.github.io/NuGridPy) is the NuGrid Python package providing tools to access and analyse (e.g. plot) various output formats (including the [NuGrid hdf5-based se format](https://github.com/NuGrid/NuSE)) from NuGrid codes (mppnp and ppn) and from the MESA stellar evolution code using the [NuGrid mesa_h5 MESA extension](https://github.com/NuGrid/mesa_h5). In principle the NuGridPy package can be used on any stellar evolution code output if the [NuGrid se libraries](https://github.com/NuGrid/NuSE) are used for output.<br>
The mesa.py module will work with MESA ASCII output in the `LOGS` directory.


## Using NuGridPy

The NuGridPy modules were written with an interactive work mode in mind, either

- taking advantage of the interactive ipython session, or
- inside a jupyter ipython notebooks. Once your session starts import modules, such as `mesa`, `nugridse` or `ppn` (depending on which type of data you are working with) from the NuGridPy package, for example:

```
    from nugridpy import mesa as ms
```

### Example session
* [Star explore notebook on Jupyter nbviewer](https://nbviewer.jupyter.org/github/NuGrid/wendi-examples/blob/master/Stellar%20evolution%20and%20nucleosynthesis%20data/Star_explore.ipynb)
* [More examples on Jupyter nbviewer](https://nbviewer.jupyter.org/github/NuGrid/wendi-examples/tree/master/Stellar%20evolution%20and%20nucleosynthesis%20data/Examples)

A typical example session in a jupyter notebook that can be performed at the [Web-Exploration of NuGrid Data Interactive (WENDI)](https://wendi.nugridstars.org) server would look like this:

1. Go to https://wendi.nugridstars.org and sign-in with your github ID (sessions will be culled at regular intervals > a few hours, if you want to use this service beyond this trial period send a message to fherwig at uvic.ca)
2. Start a Python 3 ipython notebooks
3. Load NuGridPy packages and initialise data source:
```
%pylab
# loading packages
from nugridpy import nugridse as nuse
from nugridpy import mesa
#setting data path for mesa and nuse
data_dir='/data/nugrid_vos'
# data_dir='/data/nugrid_apod2/' # alternative data store
# do ! ls /data/nug* to check for other alternative data stores
mesa.set_nugrid_path(data_dir)
nuse.set_nugrid_path(data_dir)
```
4. Creating see and ppd instances
```
# see: Stellar Evolution and Explosion data
# ppd: Post-Processing Data
m2z02_ppd=nuse.se(mass=2,Z=0.02)
m2z02_see=mesa.history_data(mass=2,Z=0.02)
```
5. Plot Hertzsprung-Russel diagram or Kippenhahn diagram
```
m2z02_see.hrd_new()
```
```
m2z02_see.kip_cont()
```
6. Inquire doc string, plot abundance profiles from ppd data_dir
```
m2z02_see.plot?
```
```
figure(11)
m2z02_ppd.plot('mass','Ba-138',fname=33500,logy=True,shape='-',\
               linewidth=2,limits=[0.5882, 0.5889,-7.8, -3.2])
```
```
species=['H-1','C-12','C-13','N-14','Fe-56','Sr-86','Ba-138','Pb-206']
ifig=121;close(ifig);figure(ifig)
m2z02_ppd.abu_profile(isos=species, ifig=ifig, fname=45500, logy=True, colourblind=True)
ylim(-9,0)
xlim(0.603,0.6033)
title("Formation of the $^\mathsf{13}\mathsf{C}$ pocket: the partial H-$^\mathsf{12}\mathsf{C}$ zone")
```

## Documentation

Each module, class, function has (or should have!) reasonable doc strings in the modules. Read the docstring: `help(ms)`, `m2z02_see.plot?`

The docstrings are also available on the [Documentation web page](https://nugrid.github.io/NuGridPy/documentation.html).

If you have made tested and debugged improvements we are happy to know about them. Make a pull request on github. Such improvements include the documentation.

The tools provided here are useful to us, but of course there are still many things that need attention and improvement. Please add submit an _issue_ on the github repo for ideas of improvements and to report bugs. Let us know if you want to help with these. Pull requests and new issues are most welcome!


## Installation

There are several ways you can install NuGridPy.

### PyPI
Major release from PyPI:
```
pip install nugridpy
```

### Release from github:
Sometimes you want to install a specific release. Go to the [NuGridPy Release page](https://github.com/NuGrid/NuGridPy/releases) and determine the tag of the release you want. If the tag is `v0.7.2` install that release with pip using the following (you could choose something else for the egg name):
```
pip install -e git://github.com/NuGrid/NuGridPy.git@v0.7.4#egg=nugridpy
```

If you just want to install whatever the latest commit is on github using github you can do:
```
pip install git+https://github.com/NuGrid/NuGridPy.git
```

### Clone and PYTHONPATH
Especially for developing NuGridPy you may want to use pip but have more control of where the installation goes, to change the repo and commit back. In that case you can clone this repo, e.g.:
```
cd ; mkdir src; cd src
git clone https://github.com/NuGrid/NuGridPy.git
```
and point the `PYTHONPATH` variable to the NuGridPy repo directory.

Inside a jupyter notebook you can set the path the following way:
```
import sys
sys.path.append('/home/user/src/NuGridPy')
```


### Required packages

All modules should work with the python distribution recommended [NuGridDoc python](https://github.com/NuGrid/NuGridDoc/blob/master/Resources/Python.md) distribution, with one additional package, the _future_ package that needs to be installed additionally.

NuGridPy has the following python dependencies:
`numpy scipy matplotlib h5py xlrd future`

For additional details on required packages, dependencies and manual installation please consult the Wiki.

