PoSB-strand-displacement
========================

MIT Principles of Synthetic Biology Term Project: Modeling in vivo strand displacement cascades to identify regimes of digital logic robustness to rates of RNA production and degradation.

Files
-----

analyze_system.py -- a few examples showing how to cycle through inputs to AND/OR/NOT and record the output levels after a certain amount of time in a convenient table.

Gates.py -- classes defining the AND, OR, NOT, INP, REP components

helpers.py -- a few helper functions used here and there

main.py -- the main System class is defined here for now. Has examples of how to define AND/OR/NOT systems

Strands.py -- classes defining the strand-based components (SignalStrand, GateStrand, GO/Gate-Output complex, DynamicGate,  Th/Threshold, Reporter

sys_AND.py -- the output code of the AND system, allows one to analyze an AND gate

sys_NOT.py -- similar for the NOT gate

sys_OR.py -- similar for the OR gate
