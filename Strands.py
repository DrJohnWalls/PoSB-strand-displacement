from helpers import *

"""							Strand-based components							"""

class StrandComponent():
	def __hash__(self):
		return hash(self.rxn_name)
	def __repr__(self):
		return self.__str__()

class SignalStrand(StrandComponent):
	"""
	Strands like
	3' left domain - T - right domain 5'
	"""
	
	fuel_domain = "f"
	
	def __init__(self, left_domain, right_domain, rel_amount = 1.0):
		self.left_domain = str(left_domain)
		self.right_domain = str(right_domain)
		self.rel_amount = rel_amount
		self.rxn_name = "%s_T_%s" % \
			(format_domain_rxn(self.left_domain), format_domain_rxn(self.right_domain))
			# 3' -> 5' description
	
	def __str__(self):
		return "%s\tT\t%s" % \
			(format_domain(self.left_domain), format_domain(self.right_domain))
		# this is 3' -> 5'

class GateStrand(StrandComponent):
	"""
	Strands like
	5' T* - main_domain* - T* 3'
	"""
	
	def __init__(self, main_domain, rel_amount = 1.0):
		self.main_domain = str(main_domain)
		self.rel_amount = rel_amount
		self.right_partners = []
		self.left_partners = []
		self.rxn_name = "Tc_%s_Tc" % format_domain_rxn(self.main_domain)
			# 5' -> 3' description
	
	def __str__(self):
		return "T*\t%s*\tT*" % format_domain(self.main_domain)
		# this is 5' -> 3'


class GO(StrandComponent):
	"""
	Complexes like
	3'        |-------Si--------|--T---|-------Sk--------| 5'
	5' |--T*--|-------Si*-------|--T*--|                   3'
	"""
	
	def __init__(self, main_domain, out_domain, rel_amount = 1.0):
		self.main_domain = main_domain
		self.out_domain = out_domain
		self.gate = GateStrand(main_domain)
		self.out_signal = SignalStrand(main_domain, out_domain)
		self.rel_amount = rel_amount
		self.rxn_name = rxn_name_components(self.out_signal, self.gate)
	
	def __str__(self):
		return "\t%s\n%s" % \
			(str(self.out_signal), str(self.gate))

class DynamicGate(StrandComponent):
	"""
	Complexes like
	3'        |-------Si--------|-------Sj--------|        5'
	5' |--T*--|-------Si*-------|-------Sj*-------|--T*--| 3'
	"""
	def __init__(self, left_domain, right_domain, rel_amount = 1.0):
		self.left_domain = left_domain
		self.right_domain = right_domain
		self.rel_amount = rel_amount
		self.right_partners = []
		self.left_partners = []
		self.rxn_name = "%s_%s-Tc_%sc_%sc_Tc" % \
			(format_domain_rxn(self.left_domain), format_domain_rxn(self.right_domain),
			format_domain_rxn(self.left_domain), format_domain_rxn(self.right_domain))
	
	def __str__(self):
		return "\t%s\t%s\nT*\t%s\t%s\tT*" % \
			(format_domain(self.left_domain), format_domain(self.right_domain),
			format_domain(self.left_domain), format_domain(self.right_domain))

class Th(StrandComponent):
	
	def __init__(self, in_toehold_domain, in_domain, rel_amount = 1.0):
		self.in_toehold_domain = in_toehold_domain
		self.in_domain = in_domain
		self.rel_amount = rel_amount
		self.rxn_name = "%s-%s_Tc_%sc" % \
			(format_domain_rxn(self.in_domain),
			format_domain_rxn(self.in_toehold_domain).lower(),
			format_domain_rxn(self.in_domain))
	
	def rxn_products(self, in_signal):
		return [ "%s-%s_Tc_%sc" % \
				(in_signal.rxn_name,
				format_domain_rxn(self.in_toehold_domain).lower(),
				format_domain_rxn(self.in_domain)),
			format_domain_rxn(self.in_domain) ]
	
	def __str__(self):
		return "\t\t%s\n%s\tT*\t%s*" % \
			(format_domain(self.in_domain),
			format_domain(self.in_toehold_domain).lower(),
			format_domain(self.in_domain))

class Reporter(StrandComponent):
	def __init__(self, main_domain, rel_amount = 1.0):
		self.main_domain = main_domain
		self.rel_amount = 1.0
		self.rxn_name = "%sQ-Tc_%scF" % \
			(format_domain_rxn(self.main_domain),
			format_domain_rxn(self.main_domain))
	
	def rxn_products(self, in_signal):
		return [ "%s-Tc_%scF" % \
				(in_signal.rxn_name,
				format_domain_rxn(self.main_domain)),
			format_domain_rxn(self.main_domain + "Q") ]
	
	def __str__(self):
		return "\t%s\tQ\nTc_%sc\tF" % \
			(format_domain_rxn(self.main_domain),
			format_domain_rxn(self.main_domain))
