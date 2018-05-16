import os
import sys
from argparse import ArgumentParser

from gooey import Gooey
from gooey import GooeyParser

#import pandas
#import wx

#from gooey_2 import main2



@Gooey(show_sidebar=True)
def main():

	#print(dir(GooeyParser))

	parser = GooeyParser(description='''Test''')

	parser.add_argument('--foo', action='store_true', help='foo help')
	subparsers = parser.add_subparsers(help='types of A')

	a_parser = subparsers.add_parser("A")
	#b_parser = subparsers.add_parser("B")

	# parser = GooeyParser(
	# 	description='''Test''')
	a_parser.add_argument(
						'-i', '--input',
						required=False,
						dest = 'csv',
						help='Test_file',
						widget='FileChooser',
						)


	args = parser.parse_args()

	import wx
	app = wx.App()
	app.Destroy()

main()

# @Gooey(show_sidebar=False)
# def main2():

# 	parser2 = GooeyParser(description='''Test2''')

# 	parser2.add_argument('-test', action='store_true', help='foo help')

# 	args2 = parser2.parse_args()

# 	print(args2)

# 	# for d in dir(): print(d)
# 	# for l in locals(): print(l)
# 	# for g in globals(): print(g)

# 	# topFrame = wx.GetApp()#.GetTopWindow() 
# 	# print(topFrame)

# main()
# main2()
#exit()



# 	p = Process(target=main2)
# 	p.start()

# main()

# 	exit()
# 	print(args2)

# main()
#main2()
#	args = parser.parse_args()

	# b_parser.add_argument(
	# 					'-o', '--output',
	# 					default = args.csv,
	# 					required=False,
	# 					dest = 'out',
	# 					help='Test_file',
	# 					widget='FileChooser',
	# 					)

	# args = parser.parse_args()

	# args = parser.parse_args()

	# if args.csv:
	# 	print(args.csv)
	# b_parser.add_argument(
	# 			'-o', '--output',
	# 			default = 'args.csv',
	# 			required=False,
	# 			dest = 'out',
	# 			help='Test_file',
	# 			choices = ['Choice 1', 'Choice 2', 'Choice 3'],
	# 			nargs='+',
	# 			widget='Listbox',
	# 			)

	# args = parser.parse_args()
	# args = parser.parse_args()

	# required = parser.add_argument_group('Required arguments')

	# required.add_argument(
	# 					'-i', '--input',
	# 					required=False,
	# 					dest = 'csv',
	# 					help='Test_file',
	# 					#widget='FileChooser',
	# 					)

	# optional = parser.add_argument_group('Optional arguments')

	# optional.add_argument(
	# 					'-o', '--output',
	# 					required=False,
	# 					dest = 'out',
	# 					help='Test_file',
	# 					choices = ['Choice 1', 'Choice 2', 'Choice 3'],
	# 					nargs='+',
	# 					#widget='Listbox',
	# 					)

	# args = parser.parse_args()

	# parser = ArgumentParser(
	# 	description='''Test''')

	# required = parser.add_argument_group('Required arguments')
	# required.add_argument(
	# 				'-t', '--test',
	# 				required=False,
	# 				dest = 'test',
	# 				help='Test_file',
	# 				#widget='FileChooser',
	# 				)

	# args = parser.parse_args()

	# print(args.csv)
	#print(args.out)
	#print(args.test)


	# required2 = parser.add_argument_group('Required arguments')

	# required2.add_argument(
	# 				'-a',
	# 				required=False,
	# 				dest = 'csv',
	# 				help='Test_file',
	# 				#widget='FileChooser',
	# 				)

	# args2 = parser2.parse_args()

	# print(args2.csv)


# if __name__ == "__main__":
# 	main()