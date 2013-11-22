# TODO: put in git

from helpers import *
from Strands import *
from Gates import *

from sys import argv
print argv[0]

"""								Main class									"""

class System():
	def __init__(self, name, rates = None, is_dynamic = True, one_x = 100):
		self.internal_domain_num = 10 - 1
		self.name = name
		
		if rates == None:
			"""
			self.rates = {			# 20C, 1 unit of concentration (C) = 1 nM
				"k_s": 5e4 * 1e-9,	# see-sawing							C^-1 s^-1
				"k_f" : 2e6 * 1e-9,	# thresholding							C^-1 s^-1
				"k_rs" : 0.5,		# threshord universal toehold binding	s^-1
				"k_rf" : 10,		# regular universal toehold binding		s^-1
				"k_l" : 1 * 1e-9	# leak									C^-1 s^-1
			}
			"""
			self.rates = {			# 25C, 1 unit of concentration (C) = 1 nM
				"k_s": 5e4 * 1e-9,	# see-sawing							C^-1 s^-1
				"k_f": 2e6 * 1e-9,	# thresholding							C^-1 s^-1
				"k_rs": 1.3,		# threshord universal toehold binding	s^-1
				"k_rf": 26,		# regular universal toehold binding		s^-1
				"k_l": 10 * 1e-9,	# leak									C^-1 s^-1
				
				"transcription": 0.1,
				"degradation": 0.001,
				"dsdegradation": 0.001
			}
		else:
			self.rates = rates
		
		self.components = []
		self.inputs = []
		self.reporters = []
		self.is_dynamic = is_dynamic
		self.one_x = one_x
	
	def new_domain(self):
		self.internal_domain_num += 1
		return "s%d" % self.internal_domain_num
	
	def add_component(self, component_type, *args):
		new_component = component_type(self, *args)
		self.components.append(new_component)
		if component_type == INP:
			self.inputs.append(new_component)
		if component_type == REP:
			self.reporters.append(new_component)
		return new_component
	
	def get_reactions(self):
		# TODO: not safe to call this more than once! the part where it sets partners is problematic
		rxns = []
		all_species = set()
		SignalStrands, GOs, Ths, DynamicGates = set(), set(), set(), set()
		TranscribedSignals = set()
		GateStrands = set()
		Reporters = set()
		# TODO: perhaps getting these things should be functions of the original classes?
		for c in self.components:
			if c.gate_type == "AND" or c.gate_type == "OR":
				SignalStrands.update(set([ c.summer.out_signal, c.amplifier.out_signal, c.fuel ]))
				TranscribedSignals.add(c.fuel)
				GOs.update(set([ c.summer, c.amplifier ]))
				GateStrands.update(set([ c.summer.gate, c.amplifier.gate ]))
				Ths.add(c.threshold)
			elif c.gate_type == "NOT":
				DynamicGates.add(c.dynamic_gate)
				SignalStrands.update(set([ c.carrier, c.fuel, c.buffer.out_signal ]))
				TranscribedSignals.update(set([ c.carrier, c.fuel ]))
				GOs.add(c.buffer)
				GateStrands.add(c.buffer.gate)
			#elif c.gate_type == "BUF": # TODO: remove, a buffer is like an OR gate with a single input
			#	SignalStrands.update(set([ c.fuel, c.summer.out_signal ]))
			#	TranscribedSignals.add(c.fuel)
			#	GOs.add(c.summer)
			#	GateStrands.add(c.summer.gate)
			#	Ths.add(c.threshold)
			#
			elif c.gate_type == "INP":
				SignalStrands.add(c.signal)
				TranscribedSignals.add(c.signal)
			elif c.gate_type == "REP":
				Reporters.add(c.reporter)
			else:
				raise RuntimeError("Unknown logic gate type %s" % c.gate_type)
		
		# seesaw reactions
		# first find all possible partners a gate strand can have
		for s in SignalStrands:
			for g in GateStrands:
				if s.right_domain == g.main_domain:
					g.left_partners.append(s) # this strand can be attached on the left
				if s.left_domain == g.main_domain: # not elif!
					g.right_partners.append(s) # this strand can be attached on the right
			for dg in DynamicGates:
				if s.right_domain == dg.left_domain:
					dg.left_partners.append(s)
				if s.left_domain == dg.right_domain:
					dg.right_partners.append(s)
		# then get the reactions
		for g in GateStrands:
			for i in g.left_partners: # inputs attach from left
				for o in g.right_partners: # outputs attach from right
					substrates = [ i.rxn_name, rxn_name_components(o, g) ]
					products = [ rxn_name_components(i, g), o.rxn_name ]
					rxns.append(Reaction(substrates, products,
						"k_s", "see-saw"))
					rxns.append(Reaction(products, substrates,
						"k_s", "see-saw (rev)"))
					# because we do the reverse reaction here as well, no need
					# to loop through right parterns and then left partners
					all_species.update(set(substrates))
					all_species.update(set(products))
		for dg in DynamicGates:
			# TODO: doesn't _really_ support multiple inputs here? would have to model like GO, where only the bottom strand is looked at.
			for l in dg.left_partners: # reaction from the left first
				substrates = [ l.rxn_name, dg.rxn_name ]
				products = [ "%s-%s" % (l.rxn_name, dg.rxn_name) ]
				rxns.append(Reaction(substrates, products,
					"k_s", "see-saw"))
				rxns.append(Reaction(products, substrates,
					"k_s", "see-saw (rev)"))
				# because we do the reverse reaction here as well, no need
				# to loop through right parterns and then left partners
				all_species.update(set(substrates))
				all_species.update(set(products))
				for r in dg.right_partners: # while something is bound on the left, something can attach on the right as well
					substrates = [ "%s-%s" % (l.rxn_name, dg.rxn_name), r.rxn_name ]
					products = [ "%s-%s-%s" % (l.rxn_name, dg.rxn_name, r.rxn_name) ]
					rxns.append(Reaction(substrates, products,
						"k_s", "see-saw (irreversible cooperative hybridization)"))
						# TODO: output strand?
					all_species.update(set(substrates))
					all_species.update(set(products))
			for r in dg.right_partners: # reaction from the right first
				substrates = [ r.rxn_name, dg.rxn_name ]
				products = [ "%s-%s" % (dg.rxn_name, r.rxn_name) ]
				rxns.append(Reaction(substrates, products,
					"k_s", "see-saw"))
				rxns.append(Reaction(products, substrates,
					"k_s", "see-saw (rev)"))
				# because we do the reverse reaction here as well, no need
				# to loop through right parterns and then left partners
				all_species.update(set(substrates))
				all_species.update(set(products))
				for l in dg.left_partners: # while something is bound on the right, something can attach on the left as well
					substrates = [ "%s-%s" % (dg.rxn_name, r.rxn_name), l.rxn_name ]
					products = [ "%s-%s-%s" % (l.rxn_name, dg.rxn_name, r.rxn_name) ]
					rxns.append(Reaction(substrates, products,
						"k_s", "see-saw (irreversible cooperative hybridization)"))
					all_species.update(set(substrates))
					all_species.update(set(products))
		
		# production reactions
		## first, take into account all input-output counts by updating the
		## relevant rel_amounts parameters
		for c in self.components:
			c.set_rel_amounts()
		## find all species to be produced
		transcribed_species_ratios = set()
		for c in TranscribedSignals | GOs | Ths | DynamicGates:
			# note: reporters aren't transcribed!
			transcribed_species_ratios.add((c.rxn_name, c.rel_amount))
			all_species.add(c.rxn_name)
		if self.is_dynamic:
			for t, ra in transcribed_species_ratios:
				rxns.append(Reaction([ "" ], [ t ], "%s * transcription" % ra, comment = "production"))
			
		
		# thresholding reactions
		for s in SignalStrands:
			for th in Ths:
				if th.in_domain == s.right_domain and th.in_toehold_domain == s.left_domain:
					substrates = [ s.rxn_name, th.rxn_name ]
					products = th.rxn_products(s)
					rxns.append(Reaction(substrates, products,
						"k_f", comment = "thresholding"))
					all_species.update(set(substrates))
					all_species.update(set(products))
		
		non_degrading = set()
		# reporting reactions
		for r in Reporters:
			for s in SignalStrands:
				if s.right_domain == r.main_domain:
					substrates = [ s.rxn_name, r.rxn_name ]
					products = r.rxn_products(s)
					rxns.append(Reaction(substrates, products,
						"k_s", comment = "reporting"))
					all_species.update(set(substrates))
					all_species.update(set(products))
					non_degrading.update(set(products))
					non_degrading.add(r.rxn_name)
		
		# side reactions: universal toehold binding reactions
		# TODO: toehold leak even for correct partners?
		# TODO: universal toehold leak with reporters
		# TODO: universal toehold binding reactions for DynamicGates
		# TODO: "bulk" leak reaction option
		for s in SignalStrands:
			for g in GateStrands:
				for a in g.right_partners:
					substrates = [ s.rxn_name, rxn_name_components(a, g) ]
					products = [ "%s-%s" % (s.rxn_name, rxn_name_components(a, g)) ]
					rxns.append(Reaction(substrates, products,
						"k_f", "universal toehold leak: from left"))
					rxns.append(Reaction(products, substrates,
						"k_rf", "universal toehold leak: from left (rev)"))
					all_species.update(set(substrates))
					all_species.update(set(products))
				for a in g.left_partners:
					substrates = [ s.rxn_name, rxn_name_components(a, g) ]
					products = [ "%s-%s" % (rxn_name_components(a, g), s.rxn_name) ]
					rxns.append(Reaction(substrates, products,
						"k_f", "universal toehold leak: from right"))
					rxns.append(Reaction(products, substrates,
						"k_rf", "universal toehold leak: from right (rev)"))
					all_species.update(set(substrates))
					all_species.update(set(products))
			for th in Ths:
				substrates = [ s.rxn_name, th.rxn_name ]
				products = [ "%s-%s" % (s.rxn_name, th.rxn_name) ]
				rxns.append(Reaction(substrates, products,
					"k_f", "universal toehold leak: threshold" ))
				rxns.append(Reaction(products, substrates,
					"k_rs", "universal toehold leak: threshold (rev)" ))
				all_species.update(set(substrates))
				all_species.update(set(products))
		# side reactions: leak reactions
		# TODO: for dynamic gates, thresholds?
		for g in GateStrands:
			for r1 in g.right_partners:
				for r2 in g.right_partners:
					substrates = [ r1.rxn_name, rxn_name_components(r2, g) ]
					products = [ r2.rxn_name, rxn_name_components(r1, g) ]
					rxns.append(Reaction(substrates, products,
						"k_l", "right-leak"))
					rxns.append(Reaction(products, substrates,
						"k_l", "right-leak (rev)"))
					all_species.update(set(substrates))
					all_species.update(set(products))
			for l1 in g.left_partners:
				for l2 in g.left_partners:
					substrates = [ l1.rxn_name, rxn_name_components(l2, g) ]
					products = [ l2.rxn_name, rxn_name_components(l1, g) ]
					rxns.append(Reaction(substrates, products,
						"k_l", "left-leak"))
					rxns.append(Reaction(products, substrates,
						"k_l", "left-leak (rev)"))
					all_species.update(set(substrates))
					all_species.update(set(products))
		
		# degradation reactions
		if self.is_dynamic:
			# TODO: should reporters be degraded? with what rate? (currently they aren't)
			for species in all_species:
				if species in non_degrading: # reporters and associated complexes won't get degraded
					continue
				if "-" in species: # a complex
					rxns.append(Reaction([ species ], [ "" ], "dsdegradation", comment = "ds degradation"))
				else: # single strand species
					rxns.append(Reaction([ species ], [ "" ], "degradation", comment = "degradation"))
		return rxns, all_species, transcribed_species_ratios
	
	def LBS(self, initial_concentrations):
		rxns, species, transcribed_species_ratios = self.get_reactions()
		r = ""
		r += "directive sample 43200.0 1000\n"
		r += "directive concentration nM\n"
		r += "\n"
		
		for rate_constant in self.rates:
			r += "rate %s = %g;\n" % (rate_constant, self.rates[rate_constant])
		r += "\n"
		
		if not self.is_dynamic:
			for species, ratio in transcribed_species_ratios:
				r += "init %s %g | \n" % (species, ratio * self.one_x)
		
		else:
			for inp in self.inputs:
				r += "rate %s = ; (* relative amount of input %s *)\n" % (inp.signal.rel_amount, inp.signal.rxn_name)
		
		r += "\n"
		
		for species in initial_concentrations:
			r += "init %s %g |\n" % (species, initial_concentrations[species])
		r += "\n"
		
		r += "\t|\n".join([ str(r) for r in rxns ])
		return r
	
	def python(self, initial_concentrations):
		rxns, species, transcribed_species_ratios = self.get_reactions()
		all_species = sorted(list(species))
		
		num_species = len(all_species)
		
		ind = dict( [ (all_species[i], i) for i in xrange(num_species) ] )
		
		r = ""
		# TODO: function for getting endpoint only?
		r += "from scipy import integrate\n"
		r += "from numpy import linspace\n\n"
		r += "def analyze_system(max_t, num_samples, " + ", ".join([ inp.signal.rel_amount for inp in self.inputs ]) + ", get_final_states):\n"
		r += "\t\"\"\" Returns the final state of interesting species for the system \"%s\" \"\"\"\n" % self.name
		r += "\tall_species = %s\n" % all_species
		r += "\tnum_species = len(all_species)\n\t\n"
		r += "\tinitial = [ 0.0 ] * num_species # initial conditions\n"
		for init_species in initial_concentrations:
			r += "\tinitial[%d] = %g # %s\n" % (ind[init_species], initial_concentrations[init_species], init_species)
		if self.is_dynamic:
			for inp in self.inputs:
				r += "\tinitial[%d] = %s * %g # %s\n" % (ind[inp.signal.rxn_name], inp.signal.rel_amount, self.one_x, inp.signal.rxn_name)
		if not self.is_dynamic:
			for species, ratio in transcribed_species_ratios:
				r += "\tinitial[%d] = %s * %g # %s\n" % (ind[species], ratio, self.one_x, species)
		r += "\t\n"
		r += "\t# rate constants:\n"
		for rate in self.rates:
			r += "\t%s = %g\n" % (rate, self.rates[rate])
		r += "\t\n\tdef dXdY(X, t):\n"
		r += "\t\t\"\"\" Returns the gradient of the vector-valued function X \"\"\"\n"
		r += "\t\t# reaction propensities:\n"
		relevant_propensities = [ [] for i in xrange(num_species) ]
		i = 0
		for rxn in rxns:
			if rxn.substrates[0] == "":
				substrates = ""
			else:
				substrates = "".join([ " * X[%d]" % ind[substrate] for substrate in rxn.substrates ])
			r += "\t\tr_%d = %s%s # (%s)\n" % (i, rxn.fwdrate, substrates, rxn)
			for s in rxn.substrates:
				if s != "":
					relevant_propensities[ind[s]].append( "-r_%d" % i )
			for p in rxn.products:
				if p != "":
					relevant_propensities[ind[p]].append( "+r_%d" % i )
			i += 1
		r += "\t\t\n\t\t# species derivatives:\n"
		i = 0
		for specimen in all_species:
			r += "\t\tdsp_%d = %s\n" % (i, " ".join(relevant_propensities[ind[specimen]]))
			i += 1
		r += "\t\t\n\t\treturn [ %s ]\n" % ", ".join([ "dsp_%d" % i for i in xrange(num_species) ])
		r += "\t\n"
		
		r += "\tt = linspace(0, max_t, num = num_samples)\n"
		r += "\tvalues = integrate.odeint(dXdY, initial, t)\n"
		r += "\t\n"
		
		r += "\tret = []\n"
		r += "\t\n"
		r += "\t# get the final states of the species asked for\n"
		r += "\tfor specimen in get_final_states:\n"
		r += "\t\tret.append(values[-1:,all_species.index(specimen)][0])\n"
		r += "\treturn ret # return final states"
		r += "\n\n"
		r += "#print analyze_system(12 * 60 * 60, 1000, [ " + ", ".join([ format_domain_rxn(rep.main_domain) + "Q" for rep in self.reporters ]) + " ])" # TODO: handle this in a nicer (more general) way
		return r
	
	def __str__(self):
		r = "==== %s ====\n" % self.name
		return r + "\n".join([ str(c) for c in self.components ])
	def __repr__(self):
		return self.__str__()



"""								Reactions									"""
class Reaction():
	def __init__(self, substrates, products, fwdrate, comment = ""):
		self.substrates = substrates
		self.products = products
		self.fwdrate = fwdrate
		self.comment = comment
	
	def __str__(self):
		comment = ""
		if self.comment != "":
			comment = "\t(* %s *)" % self.comment
			
		fwd = "%s ->{%s} %s" % (" + ".join(self.substrates), self.fwdrate, " + ".join(self.products))
		return fwd + comment
	
	def __repr__(self):
		return self.__str__()

"""
rel_oscillator_1 = System("Relaxation oscillator, no buffer")
A = rel_oscillator_1.add_component(AND, "A", "A")
B = rel_oscillator_1.add_component(NOT, "A", "A")
print rel_oscillator_1


rel_oscillator_2 = System("Relaxation oscillator, with buffer")
A = rel_oscillator_2.add_component(AND, "A", "A")
B = rel_oscillator_2.add_component(BUF, "A", rel_oscillator_2.new_domain())
C = rel_oscillator_2.add_component(NOT, B.out_domain, "A")
print sys
"""

sys_AND = System("AND gate", is_dynamic = False)
sys_AND.add_component(INP, "a", "AND", "A_amount")
sys_AND.add_component(INP, "b", "AND", "B_amount")
sys_AND.add_component(AND, "AND", "out")
sys_AND.add_component(REP, "out")

sys_OR = System("OR gate", is_dynamic = False)
sys_OR.add_component(INP, "a", "OR", "A_amount")
sys_OR.add_component(INP, "b", "OR", "B_amount")
sys_OR.add_component(OR, "OR", "out")
sys_OR.add_component(REP, "out")

sys_NOT = System("NOT gate", is_dynamic = False)
sys_NOT.add_component(INP, "a", "NOT", "A_amount")
sys_NOT.add_component(NOT, "NOT", "out")
sys_NOT.add_component(REP, "out")

initial = {
	"SoutQ-Tc_SoutcF": 150
}
#print sys.LBS(initial)
to_file("sys_AND.py", sys_AND.python(initial))
to_file("sys_OR.py", sys_OR.python(initial))
#print sys_NOT.LBS(initial)
to_file("sys_NOT.py", sys_NOT.python(initial))
