import sys
from cx_Freeze import setup, Executable
#from distutils.core import setup

# get the gooey figs and lang
import os
import gooey

#call from sublime
sys.argv.append("build")

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

shared_libs = [
				('/usr/lib/x86_64-linux-gnu/libpng12.so.0','lib/wx/libpng12.so.0'), #if libs are missing, just add it !
                ('/usr/lib/x86_64-linux-gnu/libSDL-1.2.so.0','lib/wx/libSDL-1.2.so.0'),
                #('/usr/lib/libBLT.2.5.so.8.6','lib/libBLT.2.5.so.8.6')  #needs to be added and tested !
                # ('/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21', 'lib/wx/libstdc++.so.6'),
                # ('/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21', 'lib/libstdc++.so.6'),
                # ('/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21', 'lib/x86_64-linux-gnu/libstdc++.so.6'),
                ('/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21', 'libstdc++.so.6'),
                # ('/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21', 'platforms/libstdc++.so.6'),
				]

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {
 "packages": ["os"], 
 'includes':[
 'numpy.core._methods',
 'numpy.lib.format',
 'scipy.sparse.csgraph._validation',
 #'wx.adv'
 #'matplotlib.backends.backend_tkagg',
 #"tkinter",
 ],
 'include_files': get_resources() + shared_libs, 
 "excludes": ["PyQt4"]
 }

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

#additional_mods = ['numpy.core._methods', 'numpy.lib.format']

setup(  name = "ssrviz",
        version = "0.1",
        description = "Subfamily specific residue (ssr) detection and visualization toolbox",
        options = {"build_exe": build_exe_options},
        url = 'http://phabi.de/',
        author='Paul Zierep',
      	author_email='Paul.Zierep@googlemail.com',
      	#executables = [Executable("SSP-viz-draw.py", base=base)])
		executables = [Executable("ssrviz_draw.py", base=base), Executable("ssrviz.py", base=base)])
        #executables = [Executable("build_csv_gui.py", base=base), Executable("plot_dpssm_gui.py", base=base)])