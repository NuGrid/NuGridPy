from setuptools import setup, find_packages
import os

def _get_version():
    """"Convenient function to get the version of this package."""

    ns = {}
    version_path = os.path.join('nugridpy','version.py')
    if not os.path.exists(version_path):
        return None
    with open(version_path) as version_file:
        exec(version_file.read(), ns)

    return ns['__version__']


setup(name='NuGridpy',
      version=_get_version(),
      description='Python tools for NuGrid',
      author='NuGrid Team',
      author_email='fherwig@uvic.ca',
      install_requires = ["numpy", "setuptools","future", "matplotlib","scipy","h5py","xlrd"],
      url='https://nugrid.github.io/NuGridPy',
      py_modules=['ppn', 'nugridse', 'mesa', 'data_plot', 'utils',
                  'ascii_table', 'h5T','grain', 'astronomy', 'constants',
                  'io', 'data'],
      # https://pypi.python.org/pypi?:action=list_classifiers : PyPI classifiers.
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console', 'Framework :: IPython',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: BSD License',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.5',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Visualization'],
      packages = find_packages(),
      license='BSD 3-clause',
      platforms=['Linux', 'OS X'],
      data_files=[('./NuGridpy', ['LICENSE', 'README.md', 'AUTHORS'])])
