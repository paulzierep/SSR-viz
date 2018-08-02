import sys
#from cx_Freeze import setup, Executable
#from distutils.core import setup

import setuptools
from setuptools import setup

# get the gooey figs and lang
import os
# import gooey

# def get_resources():
#     target_prefix = 'gooey'
#     source_dir = os.path.dirname(gooey.__file__)
#     subdirs = ['languages', 'images']
#     includes = []
#     for directory in subdirs:
#         path = os.path.join(source_dir, directory)
#         for file in os.listdir(path):
#             file_path = os.path.join(path, file)
#             relative_path = os.path.join(target_prefix, directory, file)
#             includes.append((file_path, relative_path))
#     return(includes)

# Dependencies are automatically detected, but it might need fine tuning.
# build_exe_options = {"packages": ["os"], 'includes':[
#  'numpy.core._methods',
#  'numpy.lib.format',
#  'matplotlib.backends.backend_tkagg'
#  ],
#  'include_files': get_resources(), 
#  #"excludes": ["PyQt4"]
#  }

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
#entry_points = {'console_scripts': ['ssrviz=ssrviz.ssrviz:main'],}
if sys.platform == "win32":
    base = "Win32GUI"

#additional_mods = ['numpy.core._methods', 'numpy.lib.format']

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
   name = "ssrviz",
    version = "0.1.2.10",
    description = "Subfamily specific residue (ssr) detection and visualization toolbox",
    #options = {"build_exe": build_exe_options},
    url = 'http://phabi.de/',
    author='Paul Zierep',
    author_email='Paul.Zierep@googlemail.com',
    #packages=['ssrviz'], <- this is automatically done by setuptools.find_packages() looking for the man

    long_description=long_description,
    long_description_content_type="text/markdown",

    #packages needed !
    install_requires=[
      'biopython>=1.70',
      'numpy>=1.14.3',
      'scipy>=1.1.0',
      'matplotlib>=2.2.2',
      'natsort>=5.3.0',
      'Gooey==1.0.0', 
      'pandas>=0.21.0'
    ],

    packages = setuptools.find_packages(exclude=['fx']),


    #scripts=['ssrviz_gui'],
    #this allows to call the package from the command line
    #entry_points = {} 
    entry_points = {'gui_scripts': ['ssrviz=ssrviz.ssrviz_sub:main'],}

    )

