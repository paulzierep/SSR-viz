import re

# p = re.compile('(\d{4}) fsd')
# m = p.findall('1234 fsdfsd 1236')
# print(m)
# exit()


INPUT_BIB = "bib.txt"

OUTPUT_BIB = "conv_{0}".format(INPUT_BIB)

new_text = ""

#p = re.compile("\((\d{4})\)")
#p = re.compile("\(\d{4}\) ([^.]*\.)")
p = re.compile("([^,]*)\,.*\((\d{4})\) ([^.]*\.) .* (\d*.\d*)\.")

name_index = {}

with open(INPUT_BIB) as bib_file:
	for line in bib_file.readlines():
		#print(line)
		m = p.search(line)
		if m:
			#print(m.groups())
			name = m.group(1)
			# print(m.group(0))
			# exit()

			#check the occurrence of the name
			if name in name_index:
				name_index[name] += 1
			else:
				name_index[name] = 0

			year = m.group(2)
			title = m.group(3).replace('.',',')
			pages = m.group(4).replace('–', '-')

			#todo
			name_ind = name + str(name_index[name])

			#todo
			bibitem = '\\bibitem[{0} {{\it et~al}}., {1}]{{{2}}} \n'.format(name, year, name_ind)
			cite = m.group(0) + '\n'

			bib_bioinf = bibitem + cite
			new_text += bib_bioinf


			# print(name)
			# print(name_index)
			# print(title)
			# print(pages)
			#print(name)

			#pass
		else:
			print('####')
			print('No match for this citation:')
			print(line.strip('\n'))
			print('####')
			# print(m.group(0))


print(new_text)