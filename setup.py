from setuptools import setup, find_packages

setup(
    name = "nugridpy",
    version = "v0.7.3",
    packages = find_packages(),
    install_requires = ["numpy", "setuptools","future",\
        "matplotlib","scipy","matplotlib","h5py","xlrd"],
    author = "NuGrid collaboration",
    author_email = "fherwig@uvic.ca")
