from setuptools import setup
from io import open

requirements = [
            'numpy',
            'numba',
            'scipy',
            'mpi4py',
            'h5py',
            'periodictable',
      ]

setup(name='pysingfel',
      description='Python version of singfel.',
      long_description=open('README.rst', encoding='utf8').read(),
      url='https://github.com/wang-zy/pysingfel',
      packages=['pysingfel'],
      scripts=['bin/radiationDamageMPI'],
      install_requires=requirements,
      zip_safe=False)
