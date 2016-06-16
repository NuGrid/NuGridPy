# Pull images from VOspace and compare images quantiatively

## Purpose
This testing branch to sketch the following capabilities:

* demonstrate how to pull test data from the VOspace
* demonstrate how images in png format can be compared quantitatively

## Useage
```
export PYTHONPATH=$PYTHONPATH:"/path/to/NuGridPy-parent-directory"
python3 -m NuGridPy.testing.test1
python2 -m NuGridPy.testing.test1
```

## Limitations
This has been tested on python 3.5.1, 2.7.10 (Fedora 23)

## Further steps
Please consider how this comparison technique of images can be used in the NuGridPy test suite.
