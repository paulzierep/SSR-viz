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

@Gooey
def main():


	parser = GooeyParser(
		description='''Creates the differnce pssm plot, based on the
					class_labela.csv file and the converted alignmet file''')

	parser.add_argument(
						'-i', '--input-csv',
						required=True,
						dest = 'csv',
						help='Input file, must be a csv file that with corresponding names to the alignment',
						widget='FileChooser',
						)

	parser.add_argument(
						'-cl', '--class_label',
						#required=True,
						default = 'Class',
						dest = 'cl',
						help='An alternative column with the class label in the csv can be specified',
						)

	parser.add_argument(
						'-a', '--alignment',
						required=True,
						dest = 'ali',
						help='Input file, must be a sequence alignment with corresponding names to the csv file',
						widget='FileChooser'
						)

	parser.add_argument(
						'-o', '--output',
						dest = 'output',
						default = 'dPSSM_plots',
						help='Output folder (default: dPSSM_plots)'
						)

	parser.add_argument(
						'-p', '--plot',
						dest = 'plot',
						default = 'dPSSM',
						help='Output csv file (default: dPSSM.pdf)'
						)

	parser.add_argument(
						'-d', '--delete',
						action='store_false',
						dest = 'delete',
						#default = True,
						help='Allows to overwrite the created plots'
						)

	parser.add_argument(
						'-b', '--best',
						dest = 'best',
						default = 10,
						help='Mark only the best x positions, good for visualization'
						)

	parser.add_argument('-w','--window', 
						help='Applies a window function of x positions to the plot(s)',
						dest = 'win',
						default = None,
						)

	parser.add_argument('-wt','--window_type', 
						help='Type of the window function to apply (default: mean), other are: max, min, std',
						dest = 'win_t',
						default = None,
						)

	parser.add_argument('-gi','--gap_importance', 
						help='Change the importance of the gap in the weight matrix (default 1)',
						dest = 'gi',
						default = 1,
						)

	parser.add_argument('-ma','--matrix', 
						help='''
						Use a replacement matrix (PAM, Blossom ....) to wight the replacement based on 
						'similarity' of the amino acids, see 
						http://biopython.org/DIST/docs/api/Bio.SubsMat.MatrixInfo-module.html for all possibilites''',
						dest = 'matrix',
						default = 'basic',
						)

	parser.add_argument('-tl','--top_label', 
						help='''Places a lable on top of each plot (not the heatmap), 
						only works well for limited number of labels (<= 14),
						everything else is too much anyway !!''',
						action='store_true',
						dest = 'tl',
						default = None,
						)

	parser.add_argument('-fs','--fig_size', 
						help='''Figsize different from A4 (29.7 cm, 21.0 cm)''',
						nargs = 2,
						dest = 'fs',
						default = None,
						)

	parser.add_argument('-cs','--chuncksize', 
						help='''Number of positions per plot (default: 100),
						use 'total' as argument to show the entire alignment in one plot''' ,
						#action='store_true',
						dest = 'cs',
						default = None,
						)

	parser.add_argument('-font','--font_size', 
						help='''Fontsize different from 8''',
						#action='store_true',
						dest = 'font',
						default = None,
						)

	parser.add_argument('-wr','--window_ratio', 
						help='''Ratio between the plot and the heatmap (default: 1 (plot) to 3 (heatmap))''',
						nargs = 2,
						#action='store_true',
						dest = 'wr',
						default = None,
						)

	parser.add_argument('-pr','--position_range', 
						help='''Staps in which the ticks for the alignment positin should be shown (default: 5)''',
						#action='store_true',
						dest = 'pr',
						default = None,
						)

	parser.add_argument('-jvp','--jv_plot', 
						help='''Creates an jalview annotation file from the plot''',
						action='store_true',
						dest = 'jv_plot',
						default = None,
						)

	# parser.add_argument('-jvhm','--jv_heatmap', 
	# 					help='''Creates an jalview annotation file from the heatmap''',
	# 					action='store_true',
	# 					dest = 'jv_heatmap',
	# 					default = None,
	# 					)

	parser.add_argument(
						'-no_cl_hm', '--no_class_label_heatmap',
						action='store_true',
						dest = 'no_cl_hm',
						default = False,
						help='Remove the class lebels from the heatmap, more then 20 are often too much'
						)


	parser.add_argument(
						'-no_hm_ava', '--heatmal_ava',
						action='store_true',
						dest = 'no_hm_ava',
						default = False,
						help='Remove all vs all representation in the heatmap'
						)


	parser.add_argument(
						'-no_pl_ava', '---plot_ava',
						action='store_true',
						dest = 'no_pl_ava',
						default = False,
						help='Remove all vs all representation in the plot'
						)


	parser.add_argument('-hm_all', 
						nargs='*', 
						help='''All vs all classes for the heatmap 
						(use the class labels as arguments to only show these classes) Example -hm_all class_A class_B
						'''
						, dest = 'hm_all', default = None)


	parser.add_argument('-hm_ova', nargs='*', help='One vs all classes for the heatmap (use the class labels as arguments to only show these classes) Example -hm_ova class_A',dest = 'hm_ova', default = None)
	parser.add_argument('-pl_all', nargs='*', help='All vs all classes for the plot (use the class labels as arguments to only show these classes) Example -pl_all class_A class_B', dest = 'pl_all', default = None)
	parser.add_argument('-pl_ova', nargs='*', help='One vs all classes for the plot (use the class labels as arguments to only show these classes) Example -pl_ova class_A', dest = 'pl_ova', default = None)




	args = parser.parse_args()

	# command_line = sys.argv

	# i = 0
	# new_line = []
	# for arg in command_line:
	# 	if i % 5 == 0:
	# 		new_line.append('\n')
	# 	new_line.append(arg)
	# 	i += 1

	# command_line = ' '.join(new_line)

	#execute the parsed arguments trough the pred_rule object 

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

	# print('''
	# 	#################################################################
	# 	Alignment file is being converted, output file:
	# 	{0}
	# 	#################################################################
	# 	'''.format(ALIGNMENT_PATH)
	# 	)

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
	pr.update_descri_dict(delimiter=',', name_field='Name', class_field=args.cl, exclude=[''])


	print('''
		#################################################################
		The plot
		-> {0}
		is created in the folder
		-> {1}
		#################################################################
		'''.format(args.plot, args.output)
		)

	# if args.best:
	# 	args.scip_gap = True

	# print(args.hm_ava)

	# exit()

	pr.update_prediction_rule('pssm_new', 
										alignment_index = True,
										#plot_all = True,
										prediction_name= args.output,
										plot_name = args.plot,
										visual = True,
										#scip_gap = args.scip_gap,
										#drop_panalty = 0.5,
										keep_data_folder = True,
										keep_plot = args.delete,
										no_alignment = True, #shot version for pssm_sa
										no_rule = True, #shot version for pssm_sa
										# until here basic infos
										#then plot specific infos
										get_best = int(args.best),
										no_hm_ava = args.no_hm_ava,
										no_pl_ava = args.no_pl_ava,
										hm_all = args.hm_all,
										pl_all = args.pl_all,
										hm_ova = args.hm_ova,
										pl_ova = args.pl_ova,
										gap_imp = float(args.gi),
										window = args.win,
										window_type = args.win_t,
										sub_matrix = args.matrix,
										top_label = args.tl,
										figsize = args.fs,
										chunksize = args.cs,
										fontsize = args.font,
										w_ratio = args.wr,
										tick_ratio = args.pr,
										jv_plot = args.jv_plot,
										drop_class_label = args.no_cl_hm,
										#command_line = command_line,
										#jv_heatmap = args.jv_heatmap,

										#hm_all = args.hm_all,
										)#scip_gap_window = 5, ) #plot_window = [270, 400])

if __name__ == "__main__":
	main()

#########################
# Examples
#########################

#./plot_dpssm.py -i ./PKS_AT_specificity_tests/NRPS_discri.csv -a ./PKS_AT_specificity_tests/NRPS_alignment.fasta -cl 'Class (more then 10)' -d -b 10 -p NRPS_gi_0_5_basic_DE -b 10 -w 10 -hm_all D E -pl_all D E -gi 0.5 -ma basic
#./plot_dpssm.py -i ./PKS_AT_specificity_tests/class_labels.csv -a ./PKS_AT_specificity_tests/temp_alignment.fasta -d -b 10 -p PKS_ver_02 -b 10 -w 10 -hm_ova -pl_ova -gi 0