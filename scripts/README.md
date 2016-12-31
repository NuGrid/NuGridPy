# scripts directory

This directory contains scripts that use NuGridPy and may be useful. Some of these scripts are candidates for integration into the NuGridPy modules.
Some script might be only python2 or python3 compatible. Scripts are not fully tested!

* `mesa` - plot and analyse mesa output
	- `2Dmaps_mesa` - plot advanced Kippenhahn style (R. Hirschi, Samuel Jones and Jacqueline den Hartogh)
* `set1data` - test integrety of set1 see and ppd data 
  	- `DataTests.ipynb`  - ipython notebook, open and access all (M,Z) SEE (star_log adn profile) and PPD data sets available in set1ext
	- `yields_data_tester.py` - make random plots of random (M,Z) pairs of see (MESA profiles, star_log) and ppd data; ipython script, use magic comand `%pylab` and then `ed` to edit and execute (Falk Herwig)
        - `vos_nugrdridse_tester.py` - a very simple test of data ont the CANFAR VOspace, python script: edit name of VOspace mount point, then do `$ python vos_nugridse_tester.py`
* `nugrid_set` - plot and analyze multiple mesa and mppnp runs and combine both outputs
        -  `nugrid_set.py` - various functions related to set1 extension 
