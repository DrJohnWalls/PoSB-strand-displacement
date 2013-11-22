from helpers import *
from Strands import *

"""								Logic gates									"""

class ANDORGate():
	def __init__(self, system, in_domain, out_domain):
		self.system = system
		self.in_domain = in_domain
		self.out_domain = out_domain
		
		# create summer, thresholder, amplifier
		self.summer = GO(self.in_domain, system.new_domain())
		self.threshold = Th(self.in_domain, self.summer.out_domain)
		self.amplifier = GO(self.summer.out_domain, self.out_domain)
		self.fuel = SignalStrand(self.amplifier.main_domain, SignalStrand.fuel_domain)
	
	def set_rel_amounts(self):
		n = len(self.summer.gate.left_partners) # number of inputs
		m = len(self.amplifier.gate.right_partners) # number of outputs
		
		self.summer.rel_amount = n
		self.amplifier.rel_amount = 1.0
		self.fuel.rel_amount = 2.0 * m
		
		self.set_th_rel_amount(n)
	
	def __str__(self):
		r = "%s -> %s -> %s\n" % (self.in_domain, self.gate_type, self.out_domain)
		r += "%s\n\n" % self.summer
		r += "%s\n\n" % self.threshold
		r += "%s\n\n" % self.amplifier
		r += "%s\n\n" % self.fuel
		return r
	def __repr__(self):
		return self.__str__()

class AND(ANDORGate):
	gate_type = "AND"
	
	def set_th_rel_amount(self, n):
		self.threshold.rel_amount = 1.1 * (n - 1 + 0.2)

class OR(ANDORGate):
	gate_type = "OR"
	
	def set_th_rel_amount(self, n):
		self.threshold.rel_amount = 1.1 * 0.6

class NOT():
	gate_type = "NOT"
	def __init__(self, system, in_domain, out_domain):
		self.system = system
		self.in_domain = in_domain
		self.out_domain = out_domain
		
		# create dynamic gate, carrier, buffer
		self.dynamic_gate = DynamicGate(self.in_domain, system.new_domain())
		self.carrier = SignalStrand(self.dynamic_gate.right_domain, system.new_domain())
		self.buffer = GO(self.carrier.right_domain, self.out_domain)
		self.fuel = SignalStrand(self.buffer.main_domain, SignalStrand.fuel_domain)
	
	def set_rel_amounts(self):
		# TODO: look at Giulio's data
		self.dynamic_gate.rel_amount = 10.0
		self.carrier.rel_amount = 0.5
		self.buffer.rel_amount = 1.0
		self.fuel.rel_amount = 2.0
	
	def __str__(self):
		r = "%s -> %s -> %s\n" % (self.in_domain, self.gate_type, self.out_domain)
		r += "%s\n\n" % self.dynamic_gate
		r += "%s\n\n" % self.carrier
		r += "%s\n\n" % self.buffer
		r += "%s\n\n" % self.fuel
		return r
	def __repr__(self):
		return self.__str__()
"""
class BUF():
	gate_type = "BUF"
	def __init__(self, system, in_domain, out_domain):
		self.system = system
		self.in_domain = in_domain
		self.out_domain = out_domain
		
		# create summer, threshold, fuel
		self.summer = GO(self.in_domain, self.out_domain)
		self.threshold = Th(self.in_domain, self.summer.out_domain)
		self.fuel = SignalStrand(self.summer.main_domain, SignalStrand.fuel_domain, SignalStrand.fuel_rel_amount)
	
	
	def __str__(self):
		r = "%s -> %s -> %s\n" % (self.in_domain, self.gate_type, self.out_domain)
		r += "%s\n\n" % self.summer
		r += "%s\n\n" % self.threshold
		r += "%s\n\n" % self.fuel
		return r
	def __repr__(self):
		return self.__str__()
"""

class INP():
	gate_type = "INP"
	def __init__(self, system, left_domain, right_domain, rel_amount):
		self.system = system
		self.signal = SignalStrand(left_domain, right_domain, rel_amount)
		self.left_domain = left_domain
		self.right_domain = right_domain
	
	def set_rel_amounts(self):
		# we already set it during initialization
		pass
	
	def __str__(self):
		r = "INP -> %s\n" % (self.signal)
		return r
	def __repr__(self):
		return self.__str__()

class REP():
	gate_type = "REP"
	def __init__(self, system, main_domain):
		self.main_domain = main_domain
		self.reporter = Reporter(main_domain)
	
	def set_rel_amounts(self):
		self.reporter.rel_amount = 1.5
	
	def __str__(self):
		r = "%s -> REP\n" % (self.main_domain)
		return r
	def __repr__(self):
		return self.__str__()
