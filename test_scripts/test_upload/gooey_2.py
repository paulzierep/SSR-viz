import os
import sys
from argparse import ArgumentParser

from gooey import Gooey
from gooey import GooeyParser

import pandas


@Gooey(show_sidebar=True)
def main2():

	parser2 = GooeyParser(description='''Test''')

	parser2.add_argument('--foo', action='store_true', help='foo help')
	subparsers2 = parser2.add_subparsers(help='types of A')

	a_parser2 = subparsers2.add_parser("B")
	#b_parser = subparsers.add_parser("B")

	# parser = GooeyParser(
	# 	description='''Test''')
	a_parser2.add_argument(
						'-t', '--input',
						required=False,
						dest = 'csv',
						help='Test_file',
						widget='FileChooser',
						)


	args2 = parser2.parse_args()
	sys.exit(main2())

main2()