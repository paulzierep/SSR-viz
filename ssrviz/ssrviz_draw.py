#!/usr/bin/env python3


import os
import argparse
#from pssm_lib.main_pred import pred_rule


#for deployment the default warnings which are
#due to biopython and notebook are switched off

import sys
import warnings

if not sys.warnoptions:
	warnings.simplefilter("ignore")

#import the prediction_rule object and give some 
#arguments via commandline

from prediction_rules.main_pred import pred_rule

from gooey import Gooey
from gooey import GooeyParser

import json
import Bio.SubsMat.MatrixInfo as MI

#SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))

def load_default_args():
	'''
	Load the def. args from a file, so that choices can be tickable
	'''
	CURRENT_PATH = os.path.dirname(sys.argv[0])
	TEMP_PATH = os.path.join(CURRENT_PATH, 'ssp_viz_temp_params.txt')
	with open(TEMP_PATH, 'r') as fp:
		data = json.load(fp) #load default args from json file

		data['matrix'] = MI.available_matrices + ['phys-chem'] #add matrices as default arguments
		#basic is removed as choice, a this is not needed with the new scoring scheme 
		return(data)

def convert_choice2list(args, unqiue_string):
	'''
	Convert 
	'''
	choice_list = []
	args_dict =  vars(args)
	for arg in args_dict:
		if unqiue_string in arg:
			if args_dict[arg]:
				cla = arg.replace(unqiue_string, '')
				choice_list.append(cla)

	return(choice_list)

@Gooey(tabbed_groups=True) #show_sidebar=True)
def main():

	#default args come form the previous window and and biopython (matrices) 
	default_args = load_default_args()

	parser = GooeyParser(
		description='''Creates the subfamily specific position plot, based on the CSV file with the
class annotation and the corresponding alignment file''')

	# parser.add_argument(
	# 					'-i', '--input-csv',
	# 					required=True,
	# 					default = default_args['csv'],
	# 					dest = 'csv',
	# 					help='Input file, must be a csv file that with corresponding names to the alignment',
	# 					widget='FileChooser',
	# 					)

	# parser.add_argument(
	# 					'-cl', '--class_label',
	# 					#required=True,
	# 					default = 'Class',
	# 					dest = 'cl',
	# 					help='An alternative column with the class label in the csv can be specified',
	# 					)

	# parser.add_argument(
	# 					'-a', '--alignment',
	# 					required=True,
	# 					dest = 'ali',
	# 					help='Input file, must be a sequence alignment with corresponding names to the csv file',
	# 					widget='FileChooser'
	# 					)

	#subs = parser.add_subparsers()

	#####################################################################################################
	#file_opts_p = subs.add_parser('File', help='File options')
	file_opts_g = parser.add_argument_group('File')

	file_opts_g.add_argument(
						'-o', '--output',
						dest = 'output',
						default = 'Data',
						help= '''Name of the output folder for the created data, this folder will
be created in the same folder as the class label CSV file
						''',
						metavar = 'Data folder',
						)

	file_opts_g.add_argument(
						'-p', '--plot',
						dest = 'plot',
						default = 'Project_01',
						help= 'Name of the created files',
						metavar = 'File name',
						)

	file_opts_g.add_argument(
						'-d', '--delete',
						action='store_true',
						dest = 'delete',
						default = False,
						help= 'Allows to overwrite the created plots',
						metavar = 'Overwrite',
						)

	file_opts_g.add_argument(
					'-v', '--verbose',
					metavar = 'Verbose Mode',
					required=False,
					dest = 'ver',
					help='Shows additional details',
					action='store_true',
					default = False,
					)

	######################################################################################################
	#algo_opts_p = subs.add_parser('Algorithm', help='Algorithm options')
	algo_opts_g = parser.add_argument_group('Algorithm')

	algo_opts_g.add_argument(
						'-ol', '--outliers',
						dest = 'outliers',
						default = 1,
						metavar = 'Outlier threshold',
						help= 'Show outlieres based on Z-score for the give threshold.'
						)

	algo_opts_g.add_argument(
						'-b', '--best',
						dest = 'best',
						default = 10,
						metavar = 'Best positions',
						help='Mark only the positions with the best X scores can be combined with the outliers'
						)

	algo_opts_g.add_argument('-gi','--gap_importance', 
						metavar = 'Gap importance',
						help='Change the importance of the gap in the weight matrix (default 0)',
						dest = 'gi',
						default = 0.1,
						)

	algo_opts_g.add_argument('-ma','--matrix', 
						help='''Use a replacement matrix (PAM, Blossom ....) to wight the replacement
based on 'similarity' of the amino acids, see http://biopython.org/DIST/docs/api/Bio.SubsMat.MatrixInfo-module.html
for an explanation.''',
						dest = 'matrix',
						widget='Dropdown',
						choices = default_args['matrix'],
						default = 'phys-chem',
						)

	algo_opts_g.add_argument('-mw','--matrix_weigth', 
						metavar = 'Weight of the replacement matrix',
						help='The matrix weight is added as an additional factor to the score',
						dest = 'mw',
						default = 0.01,
						)

	######################################################################################################
	#fig_opts_p = subs.add_parser('Figure', help='Figure options')
	fig_opts_g = parser.add_argument_group('Figure')

	fig_opts_g.add_argument('-tl','--top_label', 
					metavar = 'Remove label from top of the plot',
					help='''Remove label from top of each plot,
more then 14 are too many''',
					action='store_true',
					dest = 'tl',
					default = False,
					)

	fig_opts_g.add_argument('-fs_h','--fig_size_h', 
						metavar = 'Figsize hight',
						help='''Default: 29.7 [cm]''',
						dest = 'fs_h',
						default = 29.7,
						)

	fig_opts_g.add_argument('-fs_w','--fig_size_w', 
					metavar = 'Figsize width',
					help='''Figsize width (Default: 21.0) [cm]''',
					dest = 'fs_w',
					default = 21.0,
					)

	fig_opts_g.add_argument('-cs','--chuncksize', 
						metavar = 'Chunk size of the figure',
						help='''Number of positions per plot slide (default: 100),
use 'total' as argument to show the entire alignment in one plot''' ,
						#action='store_true',
						dest = 'cs',
						default = 100,
						)

	fig_opts_g.add_argument('-tot','--total', 
						metavar = 'Entire alignment in one plot',
						help='''Automatically adjusts the chuncksize parameter to fit the entire plot in one page''' ,
						#action='store_true',
						dest = 'tot',
						default = False,
						action='store_true',
						)


	fig_opts_g.add_argument('-font','--font_size', 
						metavar = 'Fontsize',
						help='''Fontsize different from 8''',
						#action='store_true',
						dest = 'font',
						default = 8,
						)

	fig_opts_g.add_argument('-pr','--position_range', 
						metavar = 'Tick positions',
						help='''Staps in which the ticks for the alignment position should be shown (default: 5)''',
						#action='store_true',
						dest = 'pr',
						default = 5,
						)

	fig_opts_g.add_argument('-wr_pl','--window_ratio_pl', 
						metavar = 'Ratio plot',
						help='''Ratio between the plot and the heatmap''',
						#action='store_true',
						dest = 'wr_pl',
						default = '1',
						)

	fig_opts_g.add_argument('-wr_hm','--window_ratio_hm', 
						metavar = 'Ratio heatmap',
						help='''Ratio between the plot and the heatmap''',
						#action='store_true',
						dest = 'wr_hm',
						default = '1',
						)

	######################################################################################################
	#special_opts_p = subs.add_parser('Special', help='Special options')
	special_opts_g = parser.add_argument_group('Additional output')

	special_opts_g.add_argument(
						'-jvp','--jv_plot', 
						metavar = 'Jalview annotation file',
						help='''Creates an jalview annotation file from the plot
can be loaded into Jalview as alignment annotation''',
						action='store_true',
						dest = 'jv_plot',
						default = None,
						)

	special_opts_g.add_argument(
						'-stats','--plot_stats', 
						metavar = 'Plot statistics',
						help='''Creates a csv file with plot statistics, this file 
is needed to map possible structure indices to the alignment (-> Add_pdb)''',
						action='store_true',
						dest = 'stats',
						default = None,
						)

	special_opts_g.add_argument(
						'-stats_p','--stats_percentage', 
						metavar = 'Conservation in percent to be shown in the plot statistics',
						help='''Percentage of conservation to consider for the stats (0-1)''',
						dest = 'stats_p',
						default = 0,
						)

	# parser.add_argument('-jvhm','--jv_heatmap', 
	# 					help='''Creates an jalview annotation file from the heatmap''',
	# 					action='store_true',
	# 					dest = 'jv_heatmap',
	# 					default = None,
	# 					)

	######################################################################################################
	#plot_opts_p = subs.add_parser('Plot', help='Plot options')
	plot_opts_g = parser.add_argument_group('Plot options')


	plot_opts_g.add_argument('-w','--window', 
						metavar = 'Window size',
						help='Applies a window function of x positions to the plot(s)',
						dest = 'win',
						default = 10,
						)

	plot_opts_g.add_argument('-wt','--window_type', 
						metavar = 'Window type',
						help='Type of the window function to apply on the plot',
						dest = 'win_t',
						choices = ['mean', 'max', 'min', 'std'],
						default = 'mean',
						)

	plot_opts_g.add_argument(
						'-no_pl_ava', '---plot_ava',
						metavar = 'No all vs all plot',
						action='store_true',
						dest = 'no_pl_ava',
						default = False,
						help='Remove all vs all representation in the plot'
						)

	plot_opts_all = parser.add_argument_group('One-vs-One Plot')

	classes = sorted(default_args['classes'])

	for clas in classes:
		plot_opts_all.add_argument(
			'-pl_ovo_{0}'.format(clas),
			metavar = clas,
			action='store_true',
			dest = 'pl_ovo_{0}'.format(clas))

	plot_opts_ova = parser.add_argument_group('One-vs-All Plot')

	for clas in classes:
		plot_opts_ova.add_argument(
			'-pl_ova_{0}'.format(clas),
			metavar = clas,
			action='store_true',
			dest = 'pl_ova_{0}'.format(clas))

	######################################################################################################
	#heatmap_opts_p = subs.add_parser('Heatmap', help='Heatmap options')
	heatmap_opts_g = parser.add_argument_group('Heatmap Options')

	heatmap_opts_g.add_argument(
						'-no_cl_hm', '--no_class_label_heatmap',
						metavar = 'Remove the labels from the heatmap',
						action='store_true',
						dest = 'no_cl_hm',
						default = False,
						help='Remove the class labels from the heatmap, more then 20 are often too much'
						)


	heatmap_opts_g.add_argument(
						'-no_hm_ava', '--heatmal_ava',
						metavar = 'No all vs all row in the heatmap',
						action='store_true',
						dest = 'no_hm_ava',
						default = False,
						help='Remove all vs all representation in the heatmap'
						)

	heatmap_opts_all = parser.add_argument_group('One-vs-One Heatmap')

	for clas in classes:
		heatmap_opts_all.add_argument(
			'-hm_ovo_{0}'.format(clas),
			metavar = clas,
			action='store_true',
			dest = 'hm_ovo_{0}'.format(clas))

	heatmap_opts_ova = parser.add_argument_group('One-vs-All Heatmap')

	for clas in classes:
		heatmap_opts_ova.add_argument(
			'-hm_ova_{0}'.format(clas),
			metavar = clas,
			action='store_true',
			dest = 'hm_ova_{0}'.format(clas))

	###########################
	# Parsing, default and conversion of the multiple choice fields
	###########################

	args = parser.parse_args()

	#print(args)

	args.ali = default_args['alignment']
	args.csv = default_args['csv']
	args.cl = default_args['class']
	figsize = (args.fs_h, args.fs_w)
	window_ratio = (args.wr_pl, args.wr_hm)

	pl_ovo_list = convert_choice2list(args, 'pl_ovo_')
	pl_ova_list = convert_choice2list(args, 'pl_ova_')
	hm_ovo_list = convert_choice2list(args, 'hm_ovo_')
	hm_ova_list = convert_choice2list(args, 'hm_ova_')


	# CSV_PATH = os.path.join(os.path.dirname(args.csv), args.csv)
	# ALIGNMENT_PATH = os.path.join(os.path.dirname(args.ali), args.ali)
	OUTPUT_PATH = os.path.join(os.path.dirname(args.csv))#, args.pred_rules)

	file_path_dict = 	{	
						'multi_fasta':None, 
						'pred_eval':None,
						'ref_seq':None,
						'hmmer': None,
						'pred_rules':OUTPUT_PATH,
						'alignment': args.ali, 
						'discri_csv':args.csv, 
						}

	print('''
		#################################################################
		Class label (csv) file
		-> {0}
		is being synchronised with the alignment file:
		-> {1}
		based on the csv column:
		-> {2}
		#################################################################
		'''.format(args.csv, args.ali, args.cl)
		)


	pr = pred_rule(file_path_dict, args.ali, remove_output = False) #init the object
	pr.update_descri_dict(delimiter=',', name_field='Name', class_field=args.cl, exclude=[''], reduce_DB=False, verbose = args.ver)

	print('''
		#################################################################
		The plot
		-> {0}
		is created in the folder
		-> {1}
		#################################################################
		'''.format(args.plot, args.output)
		)

	# print('########################')
	# print(args.matrix)

	#adjust the cs param.
	if args.tot:
		args.cs = 'total'

	#arguments needs to be inverted to work correctly
	args.tl = not args.tl

	#adjust the win param.
	if args.win == 0:
		args.win = None


	pr.update_prediction_rule('pssm_new', 
										alignment_index = True,
										#plot_all = True,
										prediction_name= args.output,
										plot_name = args.plot,
										visual = True,
										#scip_gap = args.scip_gap,
										#drop_panalty = 0.5,
										keep_data_folder = True,
										delete_plot = args.delete,
										no_alignment = True, #shot version for pssm_sa
										no_rule = True, #shot version for pssm_sa
										# until here basic infos
										#then plot specific infos
										get_best = int(args.best),
										get_outliers_z = float(args.outliers),
										no_hm_ava = args.no_hm_ava,
										no_pl_ava = args.no_pl_ava,

										hm_ovo = hm_ovo_list,
										pl_ovo = pl_ovo_list,
										hm_ova = hm_ova_list,
										pl_ova = pl_ova_list,

										gap_imp = float(args.gi),
										window = args.win,
										window_type = args.win_t,
										sub_matrix = args.matrix,
										top_label = args.tl,
										figsize = figsize,
										chunksize = args.cs,
										fontsize = args.font,
										w_ratio = window_ratio,
										tick_ratio = args.pr,
										jv_plot = args.jv_plot,
										stats = args.stats,
										stats_p = float(args.stats_p),
										drop_class_label = args.no_cl_hm,
										m_weigth = float(args.mw),
										verbose = args.ver,
										#command_line = command_line,
										#jv_heatmap = args.jv_heatmap,

										#hm_all = args.hm_all,
										)#scip_gap_window = 5, ) #plot_window = [270, 400])

if __name__ == "__main__":
	main()
#main()

#########################
# Examples
#########################

#./plot_dpssm.py -i ./PKS_AT_specificity_tests/NRPS_discri.csv -a ./PKS_AT_specificity_tests/NRPS_alignment.fasta -cl 'Class (more then 10)' -d -b 10 -p NRPS_gi_0_5_basic_DE -b 10 -w 10 -hm_all D E -pl_all D E -gi 0.5 -ma basic
#./plot_dpssm.py -i ./PKS_AT_specificity_tests/class_labels.csv -a ./PKS_AT_specificity_tests/temp_alignment.fasta -d -b 10 -p PKS_ver_02 -b 10 -w 10 -hm_ova -pl_ova -gi 0

#########################
#debug
#########################

# main()

# import pickle

# with open('test_pssm_dict.pickle', 'rb') as handle:
# 	pssm_dict = pickle.load(handle)

# from prediction_rules import pssm_func
# # from prediction_rules.pssm_func import all_vs_all_pssm_cons

# visual_path = ''

# sub_matrix = pssm_func.get_sub_matrix(name = 'basic', gap_importance = 0)
# df = pssm_func.all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix)
# df = pssm_func.add_all_vs_all(df) 

# pssm_func.df2pssm_visual(df = df, 
# 						path = visual_path,
# 						pssm_dict = pssm_dict,
# 						stats = True,
# 						alignment_index = True,
# 						pl_ovo = ['ema', 'mal'],
# 						pl_ova = ['mom'],
# 						get_best = 20,)	
# # print(sub_matrix)

# #sub_matrix = get_sub_matrix(name = 'blosum30', gap_importance = 0)

# #print(sub_matrix)
# exit()

#df1 = all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix, discri_keys = ['ema', 'mal'])#, 'mm', 'mom'])
# print(df1.loc[:,[310,'Class_A','Class_B']])
#print(df1.loc[:,[310,'Class_A','Class_B']])
# df1 = all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix, discri_keys = ['mal', 'ema'])#, 'mm', 'mom'])
# print(df1.loc[:,[310,'Class_A','Class_B']])

# # for k in range(10):
# #from prediction_rules.pssm_func import get_sub_matrix
# sub_matrix = get_sub_matrix(name = 'benner6', gap_importance = 1)
# #from prediction_rules.pssm_func import all_vs_all_pssm_cons
# df1 = all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix)
# print(df1.loc[:,310])

# with open('test_pssm_dict.pickle', 'rb') as handle:
# 	pssm_dict = pickle.load(handle)

# # # from prediction_rules.pssm_func import get_sub_matrix
# # # sub_matrix = get_sub_matrix(name = 'benner6', gap_importance = 0)

# # # df1 = all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix)

# from prediction_rules.pssm_func import get_sub_matrix
# from prediction_rules.pssm_func import all_vs_all_pssm_cons

# # for k in range(10):
# #from prediction_rules.pssm_func import get_sub_matrix
# sub_matrix = get_sub_matrix(name = 'benner6', gap_importance = 1)
# #from prediction_rules.pssm_func import all_vs_all_pssm_cons
# df1 = all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix)
# print(df1.loc[:,310])


#  



# print(df1)
# print(df2)
# print(df1.loc[:,310])
# print(df2.loc[:,310])

# print(df1 == df2)

# from prediction_rules.pssm_func import get_sub_matrix


# #print(sub_matrix)

# sub_matrix2 = get_sub_matrix(name = 'benner6', gap_importance = 0)
# #print(sub_matrix)

# print(sub_matrix == sub_matrix2)