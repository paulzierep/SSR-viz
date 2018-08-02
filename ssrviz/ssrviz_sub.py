#!/usr/bin/env python3

from subprocess import Popen, PIPE
import os
import sys

def main():

	CURRENT_PATH = os.path.dirname(os.path.realpath(__file__)) # <-seems to be the better option when symlinks are involved 
	SSR_VIZ_PATH = os.path.join(CURRENT_PATH, 'ssrviz.py')

	PYTHON_PATH = sys.executable
	process = Popen([PYTHON_PATH, SSR_VIZ_PATH], stdout=PIPE, stderr=PIPE)

	output, error = process.communicate()


if __name__ == "__main__":
	main()
