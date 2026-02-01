from setuptools import setup, Extension
from os import path
import numpy as np

MAJOR = 0
MINOR = 2
MICRO = 8
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

extensions = [
    Extension("pyFTracks.annealing", ["pyFTracks/annealing.pyx"]),
    Extension("pyFTracks.thermal_history", ["pyFTracks/thermal_history.pyx"]),
    Extension("pyFTracks.rci_engine", ["pyFTracks/rci_engine.pyx"]), 
]

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='pyFTracks',
    setup_requires=[
        'setuptools>=18.0',
        'numpy',
        'cython'
    ],
    version=VERSION,
    description='Fission Track Modelling and Analysis with Python',
    ext_modules=extensions,
    include_package_data=True,

    # FIX: correct way to pass NumPy include directory
    include_dirs=[np.get_include()],

    long_description=long_description,
    url='https://github.com/rbeucher/pyFTracks.git',
    author='Romain Beucher',
    author_email='romain.beucher@unimelb.edu.au',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    packages=["pyFTracks", "pyFTracks/ressources", "pyFTracks/radialplot"],
    keywords='geology thermochronology fission-tracks',
    install_requires=requirements,
)
