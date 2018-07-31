import os
import shutil

def create_folder(output_folder_path, remove = False, new_folder = True):
	'''checks if a folder exsists, if so, (remove = False) creates a new folder with increasing integer, until a a folder is found, which
	does not exsist yet, or (remove = True) replaces it with a new empty folder, returns the name of the created folder'''

	if remove == False:

		if new_folder == True:

			index = 0
			original_output_folder_path = output_folder_path

			while True:

				if os.path.isdir(output_folder_path) == False:
					os.mkdir(output_folder_path)
					print((output_folder_path + ' created'))
					break
				else:
					output_folder_path = original_output_folder_path + '_{0}'.format(index)
					index += 1

		else:
			if os.path.isdir(output_folder_path) == False:
				os.mkdir(output_folder_path)
				print((output_folder_path + ' created'))
			else:
				print((output_folder_path + ' exsists, is not replaced'))

	else:

		if os.path.isdir(output_folder_path) == True:
			shutil.rmtree(output_folder_path)
			os.mkdir(output_folder_path)
			print((output_folder_path + ' replaced'))
		else:
			os.mkdir(output_folder_path)
			print((output_folder_path + ' created'))

	return(output_folder_path)

def create_file(output_file_path, remove = False):
	'''basic file existence check and removal, return False if file is not removed and true if removed or not there''' 

	if remove == False:
		if os.path.exists(output_file_path):
			print((output_file_path + ' exsists, is not replaced'))
			return(False)
		else:
			return(True)

	else:
		if os.path.exists(output_file_path):
			os.remove(output_file_path)
			print((output_file_path + ' removed'))
		return(True)


def check_var_type(var, expected_type):
	'''checks if a variable can be transformed into a certain type, if so it is returned as this type, otherwise error'''
	if expected_type == 'int':
		try:
			var = int(var)
			return(True)
		except:
			return(False)

	if expected_type == 'str':
		try:
			var = str(var)
			return(True)
		except:
			return(False)

def check_dict_key_true(dicti, key):
	'''checks if a dict contains a certain key and if the key is True, else False'''
	if key in dicti:
		if dicti[key] == True:
			return(True)
		else:
			return(False)
	else:
		return(False)

def check_dict_key_not_none(dicti, key):
	'''checks if a dict contains a certain key and if the key is True, else False'''
	if key in dicti:
		if dicti[key] == None:
			return(False)
		else:
			return(True)
	else:
		return(False)


