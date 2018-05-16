
import os
from argparse import ArgumentParser

from gooey import Gooey
from gooey import GooeyParser

import pandas

@Gooey
def main():

	# parser = ArgumentParser(
	# 	description='''Test''')

	parser = GooeyParser(
		description='''Test''')

	required = parser.add_argument_group('Required arguments')

	required.add_argument(
						'-i', '--input',
						required=False,
						dest = 'csv',
						help='Test_file',
						widget='FileChooser',
						)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument(
						'-o', '--output',
						required=False,
						dest = 'out',
						help='Test_file',
						choices = ['Choice 1', 'Choice 2', 'Choice 3'],
						nargs='+',
						widget='Listbox',
						)

	args = parser.parse_args()

	print(args.out)
	print(args.csv)

if __name__ == "__main__":
	main()