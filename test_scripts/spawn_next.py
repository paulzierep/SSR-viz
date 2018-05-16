from gooey import Gooey
from gooey import GooeyParser

@Gooey()
def main():

	parser = GooeyParser(
		description='''Test''')

	required = parser.add_argument_group('Optional arguments')

	parser.add_argument(
						'-i', '--input',
						required=False,
						dest = 'csv',
						help='Test_file',
						widget='FileChooser',
						)


	args = parser.parse_args()

main()