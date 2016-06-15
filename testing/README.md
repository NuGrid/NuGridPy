# Pull images from VOspace and compare images quantiatively

## Purpose
This testing branch to sketch the following capabilities:

* demonstrate how to pull test data from the VOspace
* demonstrate how images in png format can be compared quantitatively

## Useage
```
export PYTHONPATH=$PYTHONPATH:"/path/to/NuGridPy"
python abu_chart.py 
python compare_image_entropy.py 
```

## Limitations
This has been tested on the Nubuntu16 VM and on the UVic scandium server using python 2.7. No VOspace authentication is needed. 

## Further steps
Please consider how this comparison technique of images can be used in the NuGridPy test suite.
