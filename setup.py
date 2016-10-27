from setuptools import setup, find_packages

setup(
    name = "nugridpy",
    version = "0.11",
    packages = find_packages(),
    install_requires = ["numpy", "setuptools","future"],
    author = "NuGrid collaboration",
    author_email = "fherwig@uvic.ca")
