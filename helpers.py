"""								Helper functions							"""

def format_domain(name):
	if len(name) == 1:
		return "S_%s" % name
	return "S_{" + name + "}"

def format_domain_rxn(name):
	return "S%s" % name


def rxn_name_components(out_signal, gate):
	return "%s-%s" % \
		(out_signal.rxn_name, gate.rxn_name)

def to_file(filename, text):
	with open(filename, 'w') as f:
		f.write(text)
