from setuptools import setup, find_packages

setup(name='NuGridpy',
      version='0.7.4', 
      description='Python tools for NuGrid',
#      long_description=open('README.rst', 'r').read(), # get discription from README
      # pandoc --from=markdown --to=rst --output=README.rst README.md
      author='NuGrid Team',
      author_email='fherwig@uvic.ca',
      install_requires = ["numpy", "setuptools","future", "matplotlib","scipy","matplotlib","h5py","xlrd"],
      url='https://nugrid.github.io/NuGridPy',
      py_modules=['ppn', 'nugridse', 'mesa', 'data_plot', 'utils',
                  'ascii_table', 'h5T','grain','astronomy'],
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
      platforms=['Linux', 'Windows', 'OS X'],
      data_files=[('./NuGridpy', ['LICENSE', 'README.md', 'AUTHORS'])])
