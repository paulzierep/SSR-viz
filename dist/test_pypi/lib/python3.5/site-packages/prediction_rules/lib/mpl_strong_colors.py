'''chosen colors for matplotlib (those which are not good for the plots are 
step by step commented !)'''

from itertools import cycle

basic_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

big_colors = ['#ffc0cb', '#708090', '#556b2f', '#ffff00', '#0000cd', '#8a2be2', '#dda0dd', '#00ffff', '#8b008b', '#006400', '#b8860b', '#dcdcdc', '#4b0082', '#ee82ee', '#cd5c5c', '#ff00ff', '#32cd32', '#6b8e23', '#008000', '#6495ed', '#808080', '#a52a2a', '#696969', '#daa520', '#1e90ff', '#ffb6c1', '#a0522d', '#00ff00', '#ffa07a', '#000000', '#ffa500', '#191970', '#7b68ee', '#228b22', '#dc143c', '#98fb98', '#800000', '#da70d6', '#9acd32', '#808000', '#adff2f', '#f0e68c', '#4169e1', '#d8bfd8', '#eee8aa', '#87cefa', '#ff69b4', '#7cfc00', '#5f9ea0', '#4682b4', '#483d8b', '#9932cc', '#7fff00', '#008b8b', '#a9a9a9', '#c0c0c0', '#bc8f8f', '#cd853f', '#fa8072', '#2f4f4f', '#00ced1', '#ff4500', '#7fffd4', '#ff8c00', '#ffdead', '#ffdab9', '#00ff7f', '#ff00ff', '#20b2aa', '#0000ff', '#c71585', '#008080', '#00fa9a', '#90ee90', '#e9967a', '#8b0000', '#8fbc8f', '#40e0d0', '#8b4513', '#000080', '#bdb76b', '#ff1493', '#ffd700', '#9400d3', '#f08080', '#778899', '#afeeee', '#ff7f50', '#48d1cc', '#faa460', '#d2b48c', '#ff6347', '#db7093', '#00008b', '#b22222', '#add8e6', '#800080', '#d2691e', '#d3d3d3', '#00ffff', '#b0e0e6', '#3cb371', '#ff0000', '#b0c4de', '#deb887', '#6a5acd', '#87ceeb', '#66cdaa', '#9370db', '#ba55d3', '#00bfff', '#2e8b57']

all_colors = cycle(basic_colors + big_colors)

def get_color_cycle():
	all_colors = cycle(basic_colors + big_colors)
	return(all_colors)

# for k in range(10):
# 	print(next(all_colors))


# new_col = []
# for color in big_colors:
# 	if color in basic_colors:
# 		pass
# 	else:
# 		new_col.append(color)

# print(len(new_col))
# print(len(big_colors))
# exit()

# print(len(big_colors_list))
# print(len(cropped_big_colors))


# print(list(big_colors_list))

# cropped_big_colors = 

# all_colors = 

# big_colors = {
# #'aliceblue':            '#F0F8FF',
# #'antiquewhite':         '#FAEBD7',
# 'aqua':                 '#00FFFF',
# 'aquamarine':           '#7FFFD4',
# #'azure':                '#F0FFFF',
# #'beige':                '#F5F5DC',
# #'bisque':               '#FFE4C4',
# 'black':                '#000000',
# #'blanchedalmond':       '#FFEBCD',
# 'blue':                 '#0000FF',
# 'blueviolet':           '#8A2BE2',
# 'brown':                '#A52A2A',
# 'burlywood':            '#DEB887',
# 'cadetblue':            '#5F9EA0',
# 'chartreuse':           '#7FFF00',
# 'chocolate':            '#D2691E',
# 'coral':                '#FF7F50',
# 'cornflowerblue':       '#6495ED',
# #'cornsilk':             '#FFF8DC',
# 'crimson':              '#DC143C',
# 'cyan':                 '#00FFFF',
# 'darkblue':             '#00008B',
# 'darkcyan':             '#008B8B',
# 'darkgoldenrod':        '#B8860B',
# 'darkgray':             '#A9A9A9',
# 'darkgreen':            '#006400',
# 'darkkhaki':            '#BDB76B',
# 'darkmagenta':          '#8B008B',
# 'darkolivegreen':       '#556B2F',
# 'darkorange':           '#FF8C00',
# 'darkorchid':           '#9932CC',
# 'darkred':              '#8B0000',
# 'darksalmon':           '#E9967A',
# 'darkseagreen':         '#8FBC8F',
# 'darkslateblue':        '#483D8B',
# 'darkslategray':        '#2F4F4F',
# 'darkturquoise':        '#00CED1',
# 'darkviolet':           '#9400D3',
# 'deeppink':             '#FF1493',
# 'deepskyblue':          '#00BFFF',
# 'dimgray':              '#696969',
# 'dodgerblue':           '#1E90FF',
# 'firebrick':            '#B22222',
# #'floralwhite':          '#FFFAF0',
# 'forestgreen':          '#228B22',
# 'fuchsia':              '#FF00FF',
# 'gainsboro':            '#DCDCDC',
# #'ghostwhite':           '#F8F8FF',
# 'gold':                 '#FFD700',
# 'goldenrod':            '#DAA520',
# 'gray':                 '#808080',
# 'green':                '#008000',
# 'greenyellow':          '#ADFF2F',
# #'honeydew':             '#F0FFF0',
# 'hotpink':              '#FF69B4',
# 'indianred':            '#CD5C5C',
# 'indigo':               '#4B0082',
# #'ivory':                '#FFFFF0',
# 'khaki':                '#F0E68C',
# #'lavender':             '#E6E6FA',
# #'lavenderblush':        '#FFF0F5',
# 'lawngreen':            '#7CFC00',
# #'lemonchiffon':         '#FFFACD',
# 'lightblue':            '#ADD8E6',
# 'lightcoral':           '#F08080',
# #'lightcyan':            '#E0FFFF',
# #'lightgoldenrodyellow': '#FAFAD2',
# 'lightgreen':           '#90EE90',
# 'lightgray':            '#D3D3D3',
# 'lightpink':            '#FFB6C1',
# 'lightsalmon':          '#FFA07A',
# 'lightseagreen':        '#20B2AA',
# 'lightskyblue':         '#87CEFA',
# 'lightslategray':       '#778899',
# 'lightsteelblue':       '#B0C4DE',
# #'lightyellow':          '#FFFFE0',
# 'lime':                 '#00FF00',
# 'limegreen':            '#32CD32',
# #'linen':                '#FAF0E6',
# 'magenta':              '#FF00FF',
# 'maroon':               '#800000',
# 'mediumaquamarine':     '#66CDAA',
# 'mediumblue':           '#0000CD',
# 'mediumorchid':         '#BA55D3',
# 'mediumpurple':         '#9370DB',
# 'mediumseagreen':       '#3CB371',
# 'mediumslateblue':      '#7B68EE',
# 'mediumspringgreen':    '#00FA9A',
# 'mediumturquoise':      '#48D1CC',
# 'mediumvioletred':      '#C71585',
# 'midnightblue':         '#191970',
# #'mintcream':            '#F5FFFA',
# #'mistyrose':            '#FFE4E1',
# #'moccasin':             '#FFE4B5',
# 'navajowhite':          '#FFDEAD',
# 'navy':                 '#000080',
# #'oldlace':              '#FDF5E6',
# 'olive':                '#808000',
# 'olivedrab':            '#6B8E23',
# 'orange':               '#FFA500',
# 'orangered':            '#FF4500',
# 'orchid':               '#DA70D6',
# 'palegoldenrod':        '#EEE8AA',
# 'palegreen':            '#98FB98',
# 'paleturquoise':        '#AFEEEE',
# 'palevioletred':        '#DB7093',
# #'papayawhip':           '#FFEFD5',
# 'peachpuff':            '#FFDAB9',
# 'peru':                 '#CD853F',
# 'pink':                 '#FFC0CB',
# 'plum':                 '#DDA0DD',
# 'powderblue':           '#B0E0E6',
# 'purple':               '#800080',
# 'red':                  '#FF0000',
# 'rosybrown':            '#BC8F8F',
# 'royalblue':            '#4169E1',
# 'saddlebrown':          '#8B4513',
# 'salmon':               '#FA8072',
# 'sandybrown':           '#FAA460',
# 'seagreen':             '#2E8B57',
# #'seashell':             '#FFF5EE',
# 'sienna':               '#A0522D',
# 'silver':               '#C0C0C0',
# 'skyblue':              '#87CEEB',
# 'slateblue':            '#6A5ACD',
# 'slategray':            '#708090',
# #'snow':                 '#FFFAFA',
# 'springgreen':          '#00FF7F',
# 'steelblue':            '#4682B4',
# 'tan':                  '#D2B48C',
# 'teal':                 '#008080',
# 'thistle':              '#D8BFD8',
# 'tomato':               '#FF6347',
# 'turquoise':            '#40E0D0',
# 'violet':               '#EE82EE',
# #'wheat':                '#F5DEB3',
# #'white':                '#FFFFFF',
# #'whitesmoke':           '#F5F5F5',
# 'yellow':               '#FFFF00',
# 'yellowgreen':          '#9ACD32'}



