#!/usr/bin/env python3

############################################################################
# This script will create a launcher for ssrviz for ubuntu/unity
############################################################################

from subprocess import Popen, PIPE
import os
import sys
from os.path import expanduser

HOME = expanduser("~")

CURRENT_PATH = os.path.dirname(os.path.realpath(__file__)) # <-seems to be the better option when symlinks are involved 
APPLICATION_PATH = os.path.join(HOME, '.local/share/applications')
DESKTOP_PATH = os.path.join(APPLICATION_PATH, 'ssrviz.desktop')

process = Popen(['which', 'ssrviz'], stdout=PIPE, stderr=PIPE)

output, error = process.communicate()

if output:
	SSRVIZ_PATH = output.decode('ascii').strip('\n')
else:
	print('ssrviz not found, is it installed ?')
	exit()

ICON_PATH = os.path.join(CURRENT_PATH,'icon.png')

desktop_file_text = '''[Desktop Entry]
Version=1.0
Name=ssrviz
Comment=ssrviz launcher
Exec={0}
Icon={1}
Terminal=false
Type=Application
Categories=Utility;Application;'''.format(SSRVIZ_PATH, ICON_PATH)

with open(DESKTOP_PATH, 'w') as desktop_file:
    desktop_file.write(desktop_file_text)


print('''{0}
	\nwritten to {1}'''.format(desktop_file_text, DESKTOP_PATH))

