#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# Some tools to help me plot using matplotlib!
import itertools as it

def line_cycle():
	return it.cycle(['-', '--', '-.', ':'])

def color_cycle():
	return it.cycle(list('bgrcmyk'))

def marker_cycle():
	return it.cycle(list('ov^<>1234sp*hH+xDd|_'))

def get_style_combos(styles_to_combine):

	styles = {'lines' : ['-', '--', '-.', ':'],
		'colors' : list('bgrcmyk'),
		'markers' : list('ov^<>1234sp*hH+xDd|_')
		}
	styles_list = [styles[st] for st in styles_to_combine]
	# print styles_list
        styles_raw = list(it.product(*styles_list))
        styles_cycle = it.cycle([''.join(x) for x in styles_raw])

	return styles_cycle

def get_style_combos_2(styles_to_combine):

	styles = {'lines' : ['-', '--', '-.', ':'],
		'colors' : list('bgrcmyk'),
		'markers' : list('ov^<>1234sp*hH+xDd|_')
		}
	styles_cycles = [it.cycle(styles[st]) for st in styles_to_combine]

	while True:
		style = ''.join([x.next() for x in styles_cycles])
		yield style

