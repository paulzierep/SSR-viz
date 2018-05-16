import sys
import os
from argparse import ArgumentParser

from subprocess import Popen, PIPE

from gooey import Gooey
from gooey import GooeyParser

import tempfile

@Gooey()
def main():

	parser = GooeyParser(
		description='''Test''')

	required = parser.add_argument_group('Optional arguments')

	parser.add_argument(
						'-i', '--input',
						required=False,
						default = 'Not there',
						dest = 'input',
						help='Test_file',
						widget='FileChooser',
						)

	args = parser.parse_args()

	########################################
	# call the next gooey as subprocess from here 
	# should work on any system
	########################################

	tmp = tempfile.NamedTemporaryFile()

	with open(tmp.name, 'w') as temp:
		temp.write(args.input)

	PYTHON_PATH = sys.executable
	process = Popen([PYTHON_PATH, 'spawn_next.py', tmp.name], stdout=PIPE, stderr=PIPE)
	output, error = process.communicate()
	print(output)
	print(error)


if __name__ == "__main__":
	main()
