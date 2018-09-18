import sys
from cx_Freeze import setup, Executable
#from distutils.core import setup

# get the gooey figs and lang
import os
import gooey

PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

def get_resources():
    target_prefix = 'gooey'
    source_dir = os.path.dirname(gooey.__file__)
    subdirs = ['languages', 'images']
    includes = []
    for directory in subdirs:
        path = os.path.join(source_dir, directory)
        for file in os.listdir(path):
            file_path = os.path.join(path, file)
            relative_path = os.path.join(target_prefix, directory, file)
            includes.append((file_path, relative_path))
    return(includes)

dlls = [
os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll'),]

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os"], 'includes':[
 'numpy.core._methods',
 'numpy.lib.format',
 'matplotlib.backends.backend_tkagg',
 'scipy.sparse.csgraph._validation',
 #'scipy.spatial.ckdtree'
 ],
 'include_files': get_resources() + dlls, 
 "excludes": ["scipy.spatial.cKDTree"]
 }

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

#additional_mods = ['numpy.core._methods', 'numpy.lib.format']

setup(  name = "ssrviz",
        version = "0.1",
        description =  "Subfamily specific residue (ssr) detection and visualization toolbox",
        options = {"build_exe": build_exe_options},
        url = 'http://www.pharmazeutische-bioinformatik.de/',
        author='Paul Zierep',
      	author_email='Paul.Zierep@googlemail.com',
      	#executables = [Executable("SSP-viz-draw.py", base=base)])
		executables = [Executable("ssrviz.py", base=base), Executable("ssrviz_draw.py", base=base)])
        #executables = [Executable("build_csv_gui.py", base=base), Executable("plot_dpssm_gui.py", base=base)])