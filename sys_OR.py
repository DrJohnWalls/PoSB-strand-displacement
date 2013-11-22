from scipy import integrate
from numpy import linspace

def analyze_system(max_t, num_samples, A_amount, B_amount, get_final_states):
	""" Returns the final state of interesting species for the system "OR gate" """
	all_species = ['SOR_T_Ss10', 'SOR_T_Ss10-SOR_T_Ss10-Tc_SOR_Tc', 'SOR_T_Ss10-Ss10-sor_Tc_Ss10c', 'SOR_T_Ss10-Ss10_T_Sf-Tc_Ss10_Tc', 'SOR_T_Ss10-Ss10_T_Sout-Tc_Ss10_Tc', 'SOR_T_Ss10-Tc_SOR_Tc', 'SOR_T_Ss10-Tc_Ss10_Tc', 'SOR_T_Ss10-Tc_Ss10_Tc-SOR_T_Ss10', 'SOR_T_Ss10-Tc_Ss10_Tc-Sa_T_SOR', 'SOR_T_Ss10-Tc_Ss10_Tc-Sb_T_SOR', 'SOR_T_Ss10-Tc_Ss10_Tc-Ss10_T_Sf', 'SOR_T_Ss10-Tc_Ss10_Tc-Ss10_T_Sout', 'SOR_T_Ss10-sor_Tc_Ss10c', 'Sa_T_SOR', 'Sa_T_SOR-SOR_T_Ss10-Tc_SOR_Tc', 'Sa_T_SOR-Ss10-sor_Tc_Ss10c', 'Sa_T_SOR-Ss10_T_Sf-Tc_Ss10_Tc', 'Sa_T_SOR-Ss10_T_Sout-Tc_Ss10_Tc', 'Sa_T_SOR-Tc_SOR_Tc', 'Sa_T_SOR-Tc_SOR_Tc-SOR_T_Ss10', 'Sa_T_SOR-Tc_SOR_Tc-Sa_T_SOR', 'Sa_T_SOR-Tc_SOR_Tc-Sb_T_SOR', 'Sa_T_SOR-Tc_SOR_Tc-Ss10_T_Sf', 'Sa_T_SOR-Tc_SOR_Tc-Ss10_T_Sout', 'Sb_T_SOR', 'Sb_T_SOR-SOR_T_Ss10-Tc_SOR_Tc', 'Sb_T_SOR-Ss10-sor_Tc_Ss10c', 'Sb_T_SOR-Ss10_T_Sf-Tc_Ss10_Tc', 'Sb_T_SOR-Ss10_T_Sout-Tc_Ss10_Tc', 'Sb_T_SOR-Tc_SOR_Tc', 'Sb_T_SOR-Tc_SOR_Tc-SOR_T_Ss10', 'Sb_T_SOR-Tc_SOR_Tc-Sa_T_SOR', 'Sb_T_SOR-Tc_SOR_Tc-Sb_T_SOR', 'Sb_T_SOR-Tc_SOR_Tc-Ss10_T_Sf', 'Sb_T_SOR-Tc_SOR_Tc-Ss10_T_Sout', 'SoutQ', 'SoutQ-Tc_SoutcF', 'Ss10', 'Ss10-sor_Tc_Ss10c', 'Ss10_T_Sf', 'Ss10_T_Sf-SOR_T_Ss10-Tc_SOR_Tc', 'Ss10_T_Sf-Ss10-sor_Tc_Ss10c', 'Ss10_T_Sf-Ss10_T_Sf-Tc_Ss10_Tc', 'Ss10_T_Sf-Ss10_T_Sout-Tc_Ss10_Tc', 'Ss10_T_Sf-Tc_Ss10_Tc', 'Ss10_T_Sout', 'Ss10_T_Sout-SOR_T_Ss10-Tc_SOR_Tc', 'Ss10_T_Sout-Ss10-sor_Tc_Ss10c', 'Ss10_T_Sout-Ss10_T_Sf-Tc_Ss10_Tc', 'Ss10_T_Sout-Ss10_T_Sout-Tc_Ss10_Tc', 'Ss10_T_Sout-Tc_SoutcF', 'Ss10_T_Sout-Tc_Ss10_Tc']
	num_species = len(all_species)
	
	initial = [ 0.0 ] * num_species # initial conditions
	initial[36] = 150 # SoutQ-Tc_SoutcF
	initial[51] = 1.0 * 100 # Ss10_T_Sout-Tc_Ss10_Tc
	initial[38] = 0.66 * 100 # Ss10-sor_Tc_Ss10c
	initial[5] = 2 * 100 # SOR_T_Ss10-Tc_SOR_Tc
	initial[24] = B_amount * 100 # Sb_T_SOR
	initial[39] = 4.0 * 100 # Ss10_T_Sf
	initial[13] = A_amount * 100 # Sa_T_SOR
	
	# rate constants:
	dsdegradation = 0.001
	k_l = 1e-08
	transcription = 0.1
	k_s = 5e-05
	k_rs = 1.3
	degradation = 0.001
	k_rf = 26
	k_f = 0.002
	
	def dXdY(X, t):
		""" Returns the gradient of the vector-valued function X """
		# reaction propensities:
		r_0 = k_s * X[24] * X[5] # (Sb_T_SOR + SOR_T_Ss10-Tc_SOR_Tc ->{k_s} Sb_T_SOR-Tc_SOR_Tc + SOR_T_Ss10	(* see-saw *))
		r_1 = k_s * X[29] * X[0] # (Sb_T_SOR-Tc_SOR_Tc + SOR_T_Ss10 ->{k_s} Sb_T_SOR + SOR_T_Ss10-Tc_SOR_Tc	(* see-saw (rev) *))
		r_2 = k_s * X[13] * X[5] # (Sa_T_SOR + SOR_T_Ss10-Tc_SOR_Tc ->{k_s} Sa_T_SOR-Tc_SOR_Tc + SOR_T_Ss10	(* see-saw *))
		r_3 = k_s * X[18] * X[0] # (Sa_T_SOR-Tc_SOR_Tc + SOR_T_Ss10 ->{k_s} Sa_T_SOR + SOR_T_Ss10-Tc_SOR_Tc	(* see-saw (rev) *))
		r_4 = k_s * X[0] * X[44] # (SOR_T_Ss10 + Ss10_T_Sf-Tc_Ss10_Tc ->{k_s} SOR_T_Ss10-Tc_Ss10_Tc + Ss10_T_Sf	(* see-saw *))
		r_5 = k_s * X[6] * X[39] # (SOR_T_Ss10-Tc_Ss10_Tc + Ss10_T_Sf ->{k_s} SOR_T_Ss10 + Ss10_T_Sf-Tc_Ss10_Tc	(* see-saw (rev) *))
		r_6 = k_s * X[0] * X[51] # (SOR_T_Ss10 + Ss10_T_Sout-Tc_Ss10_Tc ->{k_s} SOR_T_Ss10-Tc_Ss10_Tc + Ss10_T_Sout	(* see-saw *))
		r_7 = k_s * X[6] * X[45] # (SOR_T_Ss10-Tc_Ss10_Tc + Ss10_T_Sout ->{k_s} SOR_T_Ss10 + Ss10_T_Sout-Tc_Ss10_Tc	(* see-saw (rev) *))
		r_8 = k_f * X[0] * X[38] # (SOR_T_Ss10 + Ss10-sor_Tc_Ss10c ->{k_f} SOR_T_Ss10-sor_Tc_Ss10c + Ss10	(* thresholding *))
		r_9 = k_s * X[45] * X[36] # (Ss10_T_Sout + SoutQ-Tc_SoutcF ->{k_s} Ss10_T_Sout-Tc_SoutcF + SoutQ	(* reporting *))
		r_10 = k_f * X[39] * X[5] # (Ss10_T_Sf + SOR_T_Ss10-Tc_SOR_Tc ->{k_f} Ss10_T_Sf-SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left *))
		r_11 = k_rf * X[40] # (Ss10_T_Sf-SOR_T_Ss10-Tc_SOR_Tc ->{k_rf} Ss10_T_Sf + SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left (rev) *))
		r_12 = k_f * X[39] * X[29] # (Ss10_T_Sf + Sb_T_SOR-Tc_SOR_Tc ->{k_f} Sb_T_SOR-Tc_SOR_Tc-Ss10_T_Sf	(* universal toehold leak: from right *))
		r_13 = k_rf * X[33] # (Sb_T_SOR-Tc_SOR_Tc-Ss10_T_Sf ->{k_rf} Ss10_T_Sf + Sb_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_14 = k_f * X[39] * X[18] # (Ss10_T_Sf + Sa_T_SOR-Tc_SOR_Tc ->{k_f} Sa_T_SOR-Tc_SOR_Tc-Ss10_T_Sf	(* universal toehold leak: from right *))
		r_15 = k_rf * X[22] # (Sa_T_SOR-Tc_SOR_Tc-Ss10_T_Sf ->{k_rf} Ss10_T_Sf + Sa_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_16 = k_f * X[39] * X[44] # (Ss10_T_Sf + Ss10_T_Sf-Tc_Ss10_Tc ->{k_f} Ss10_T_Sf-Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_17 = k_rf * X[42] # (Ss10_T_Sf-Ss10_T_Sf-Tc_Ss10_Tc ->{k_rf} Ss10_T_Sf + Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_18 = k_f * X[39] * X[51] # (Ss10_T_Sf + Ss10_T_Sout-Tc_Ss10_Tc ->{k_f} Ss10_T_Sf-Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_19 = k_rf * X[43] # (Ss10_T_Sf-Ss10_T_Sout-Tc_Ss10_Tc ->{k_rf} Ss10_T_Sf + Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_20 = k_f * X[39] * X[6] # (Ss10_T_Sf + SOR_T_Ss10-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Tc_Ss10_Tc-Ss10_T_Sf	(* universal toehold leak: from right *))
		r_21 = k_rf * X[10] # (SOR_T_Ss10-Tc_Ss10_Tc-Ss10_T_Sf ->{k_rf} Ss10_T_Sf + SOR_T_Ss10-Tc_Ss10_Tc	(* universal toehold leak: from right (rev) *))
		r_22 = k_f * X[39] * X[38] # (Ss10_T_Sf + Ss10-sor_Tc_Ss10c ->{k_f} Ss10_T_Sf-Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold *))
		r_23 = k_rs * X[41] # (Ss10_T_Sf-Ss10-sor_Tc_Ss10c ->{k_rs} Ss10_T_Sf + Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold (rev) *))
		r_24 = k_f * X[24] * X[5] # (Sb_T_SOR + SOR_T_Ss10-Tc_SOR_Tc ->{k_f} Sb_T_SOR-SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left *))
		r_25 = k_rf * X[25] # (Sb_T_SOR-SOR_T_Ss10-Tc_SOR_Tc ->{k_rf} Sb_T_SOR + SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left (rev) *))
		r_26 = k_f * X[24] * X[29] # (Sb_T_SOR + Sb_T_SOR-Tc_SOR_Tc ->{k_f} Sb_T_SOR-Tc_SOR_Tc-Sb_T_SOR	(* universal toehold leak: from right *))
		r_27 = k_rf * X[32] # (Sb_T_SOR-Tc_SOR_Tc-Sb_T_SOR ->{k_rf} Sb_T_SOR + Sb_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_28 = k_f * X[24] * X[18] # (Sb_T_SOR + Sa_T_SOR-Tc_SOR_Tc ->{k_f} Sa_T_SOR-Tc_SOR_Tc-Sb_T_SOR	(* universal toehold leak: from right *))
		r_29 = k_rf * X[21] # (Sa_T_SOR-Tc_SOR_Tc-Sb_T_SOR ->{k_rf} Sb_T_SOR + Sa_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_30 = k_f * X[24] * X[44] # (Sb_T_SOR + Ss10_T_Sf-Tc_Ss10_Tc ->{k_f} Sb_T_SOR-Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_31 = k_rf * X[27] # (Sb_T_SOR-Ss10_T_Sf-Tc_Ss10_Tc ->{k_rf} Sb_T_SOR + Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_32 = k_f * X[24] * X[51] # (Sb_T_SOR + Ss10_T_Sout-Tc_Ss10_Tc ->{k_f} Sb_T_SOR-Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_33 = k_rf * X[28] # (Sb_T_SOR-Ss10_T_Sout-Tc_Ss10_Tc ->{k_rf} Sb_T_SOR + Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_34 = k_f * X[24] * X[6] # (Sb_T_SOR + SOR_T_Ss10-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Tc_Ss10_Tc-Sb_T_SOR	(* universal toehold leak: from right *))
		r_35 = k_rf * X[9] # (SOR_T_Ss10-Tc_Ss10_Tc-Sb_T_SOR ->{k_rf} Sb_T_SOR + SOR_T_Ss10-Tc_Ss10_Tc	(* universal toehold leak: from right (rev) *))
		r_36 = k_f * X[24] * X[38] # (Sb_T_SOR + Ss10-sor_Tc_Ss10c ->{k_f} Sb_T_SOR-Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold *))
		r_37 = k_rs * X[26] # (Sb_T_SOR-Ss10-sor_Tc_Ss10c ->{k_rs} Sb_T_SOR + Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold (rev) *))
		r_38 = k_f * X[13] * X[5] # (Sa_T_SOR + SOR_T_Ss10-Tc_SOR_Tc ->{k_f} Sa_T_SOR-SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left *))
		r_39 = k_rf * X[14] # (Sa_T_SOR-SOR_T_Ss10-Tc_SOR_Tc ->{k_rf} Sa_T_SOR + SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left (rev) *))
		r_40 = k_f * X[13] * X[29] # (Sa_T_SOR + Sb_T_SOR-Tc_SOR_Tc ->{k_f} Sb_T_SOR-Tc_SOR_Tc-Sa_T_SOR	(* universal toehold leak: from right *))
		r_41 = k_rf * X[31] # (Sb_T_SOR-Tc_SOR_Tc-Sa_T_SOR ->{k_rf} Sa_T_SOR + Sb_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_42 = k_f * X[13] * X[18] # (Sa_T_SOR + Sa_T_SOR-Tc_SOR_Tc ->{k_f} Sa_T_SOR-Tc_SOR_Tc-Sa_T_SOR	(* universal toehold leak: from right *))
		r_43 = k_rf * X[20] # (Sa_T_SOR-Tc_SOR_Tc-Sa_T_SOR ->{k_rf} Sa_T_SOR + Sa_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_44 = k_f * X[13] * X[44] # (Sa_T_SOR + Ss10_T_Sf-Tc_Ss10_Tc ->{k_f} Sa_T_SOR-Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_45 = k_rf * X[16] # (Sa_T_SOR-Ss10_T_Sf-Tc_Ss10_Tc ->{k_rf} Sa_T_SOR + Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_46 = k_f * X[13] * X[51] # (Sa_T_SOR + Ss10_T_Sout-Tc_Ss10_Tc ->{k_f} Sa_T_SOR-Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_47 = k_rf * X[17] # (Sa_T_SOR-Ss10_T_Sout-Tc_Ss10_Tc ->{k_rf} Sa_T_SOR + Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_48 = k_f * X[13] * X[6] # (Sa_T_SOR + SOR_T_Ss10-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Tc_Ss10_Tc-Sa_T_SOR	(* universal toehold leak: from right *))
		r_49 = k_rf * X[8] # (SOR_T_Ss10-Tc_Ss10_Tc-Sa_T_SOR ->{k_rf} Sa_T_SOR + SOR_T_Ss10-Tc_Ss10_Tc	(* universal toehold leak: from right (rev) *))
		r_50 = k_f * X[13] * X[38] # (Sa_T_SOR + Ss10-sor_Tc_Ss10c ->{k_f} Sa_T_SOR-Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold *))
		r_51 = k_rs * X[15] # (Sa_T_SOR-Ss10-sor_Tc_Ss10c ->{k_rs} Sa_T_SOR + Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold (rev) *))
		r_52 = k_f * X[45] * X[5] # (Ss10_T_Sout + SOR_T_Ss10-Tc_SOR_Tc ->{k_f} Ss10_T_Sout-SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left *))
		r_53 = k_rf * X[46] # (Ss10_T_Sout-SOR_T_Ss10-Tc_SOR_Tc ->{k_rf} Ss10_T_Sout + SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left (rev) *))
		r_54 = k_f * X[45] * X[29] # (Ss10_T_Sout + Sb_T_SOR-Tc_SOR_Tc ->{k_f} Sb_T_SOR-Tc_SOR_Tc-Ss10_T_Sout	(* universal toehold leak: from right *))
		r_55 = k_rf * X[34] # (Sb_T_SOR-Tc_SOR_Tc-Ss10_T_Sout ->{k_rf} Ss10_T_Sout + Sb_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_56 = k_f * X[45] * X[18] # (Ss10_T_Sout + Sa_T_SOR-Tc_SOR_Tc ->{k_f} Sa_T_SOR-Tc_SOR_Tc-Ss10_T_Sout	(* universal toehold leak: from right *))
		r_57 = k_rf * X[23] # (Sa_T_SOR-Tc_SOR_Tc-Ss10_T_Sout ->{k_rf} Ss10_T_Sout + Sa_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_58 = k_f * X[45] * X[44] # (Ss10_T_Sout + Ss10_T_Sf-Tc_Ss10_Tc ->{k_f} Ss10_T_Sout-Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_59 = k_rf * X[48] # (Ss10_T_Sout-Ss10_T_Sf-Tc_Ss10_Tc ->{k_rf} Ss10_T_Sout + Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_60 = k_f * X[45] * X[51] # (Ss10_T_Sout + Ss10_T_Sout-Tc_Ss10_Tc ->{k_f} Ss10_T_Sout-Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_61 = k_rf * X[49] # (Ss10_T_Sout-Ss10_T_Sout-Tc_Ss10_Tc ->{k_rf} Ss10_T_Sout + Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_62 = k_f * X[45] * X[6] # (Ss10_T_Sout + SOR_T_Ss10-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Tc_Ss10_Tc-Ss10_T_Sout	(* universal toehold leak: from right *))
		r_63 = k_rf * X[11] # (SOR_T_Ss10-Tc_Ss10_Tc-Ss10_T_Sout ->{k_rf} Ss10_T_Sout + SOR_T_Ss10-Tc_Ss10_Tc	(* universal toehold leak: from right (rev) *))
		r_64 = k_f * X[45] * X[38] # (Ss10_T_Sout + Ss10-sor_Tc_Ss10c ->{k_f} Ss10_T_Sout-Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold *))
		r_65 = k_rs * X[47] # (Ss10_T_Sout-Ss10-sor_Tc_Ss10c ->{k_rs} Ss10_T_Sout + Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold (rev) *))
		r_66 = k_f * X[0] * X[5] # (SOR_T_Ss10 + SOR_T_Ss10-Tc_SOR_Tc ->{k_f} SOR_T_Ss10-SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left *))
		r_67 = k_rf * X[1] # (SOR_T_Ss10-SOR_T_Ss10-Tc_SOR_Tc ->{k_rf} SOR_T_Ss10 + SOR_T_Ss10-Tc_SOR_Tc	(* universal toehold leak: from left (rev) *))
		r_68 = k_f * X[0] * X[29] # (SOR_T_Ss10 + Sb_T_SOR-Tc_SOR_Tc ->{k_f} Sb_T_SOR-Tc_SOR_Tc-SOR_T_Ss10	(* universal toehold leak: from right *))
		r_69 = k_rf * X[30] # (Sb_T_SOR-Tc_SOR_Tc-SOR_T_Ss10 ->{k_rf} SOR_T_Ss10 + Sb_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_70 = k_f * X[0] * X[18] # (SOR_T_Ss10 + Sa_T_SOR-Tc_SOR_Tc ->{k_f} Sa_T_SOR-Tc_SOR_Tc-SOR_T_Ss10	(* universal toehold leak: from right *))
		r_71 = k_rf * X[19] # (Sa_T_SOR-Tc_SOR_Tc-SOR_T_Ss10 ->{k_rf} SOR_T_Ss10 + Sa_T_SOR-Tc_SOR_Tc	(* universal toehold leak: from right (rev) *))
		r_72 = k_f * X[0] * X[44] # (SOR_T_Ss10 + Ss10_T_Sf-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_73 = k_rf * X[3] # (SOR_T_Ss10-Ss10_T_Sf-Tc_Ss10_Tc ->{k_rf} SOR_T_Ss10 + Ss10_T_Sf-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_74 = k_f * X[0] * X[51] # (SOR_T_Ss10 + Ss10_T_Sout-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left *))
		r_75 = k_rf * X[4] # (SOR_T_Ss10-Ss10_T_Sout-Tc_Ss10_Tc ->{k_rf} SOR_T_Ss10 + Ss10_T_Sout-Tc_Ss10_Tc	(* universal toehold leak: from left (rev) *))
		r_76 = k_f * X[0] * X[6] # (SOR_T_Ss10 + SOR_T_Ss10-Tc_Ss10_Tc ->{k_f} SOR_T_Ss10-Tc_Ss10_Tc-SOR_T_Ss10	(* universal toehold leak: from right *))
		r_77 = k_rf * X[7] # (SOR_T_Ss10-Tc_Ss10_Tc-SOR_T_Ss10 ->{k_rf} SOR_T_Ss10 + SOR_T_Ss10-Tc_Ss10_Tc	(* universal toehold leak: from right (rev) *))
		r_78 = k_f * X[0] * X[38] # (SOR_T_Ss10 + Ss10-sor_Tc_Ss10c ->{k_f} SOR_T_Ss10-Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold *))
		r_79 = k_rs * X[2] # (SOR_T_Ss10-Ss10-sor_Tc_Ss10c ->{k_rs} SOR_T_Ss10 + Ss10-sor_Tc_Ss10c	(* universal toehold leak: threshold (rev) *))
		r_80 = k_l * X[0] * X[5] # (SOR_T_Ss10 + SOR_T_Ss10-Tc_SOR_Tc ->{k_l} SOR_T_Ss10 + SOR_T_Ss10-Tc_SOR_Tc	(* right-leak *))
		r_81 = k_l * X[0] * X[5] # (SOR_T_Ss10 + SOR_T_Ss10-Tc_SOR_Tc ->{k_l} SOR_T_Ss10 + SOR_T_Ss10-Tc_SOR_Tc	(* right-leak (rev) *))
		r_82 = k_l * X[24] * X[29] # (Sb_T_SOR + Sb_T_SOR-Tc_SOR_Tc ->{k_l} Sb_T_SOR + Sb_T_SOR-Tc_SOR_Tc	(* left-leak *))
		r_83 = k_l * X[24] * X[29] # (Sb_T_SOR + Sb_T_SOR-Tc_SOR_Tc ->{k_l} Sb_T_SOR + Sb_T_SOR-Tc_SOR_Tc	(* left-leak (rev) *))
		r_84 = k_l * X[24] * X[18] # (Sb_T_SOR + Sa_T_SOR-Tc_SOR_Tc ->{k_l} Sa_T_SOR + Sb_T_SOR-Tc_SOR_Tc	(* left-leak *))
		r_85 = k_l * X[13] * X[29] # (Sa_T_SOR + Sb_T_SOR-Tc_SOR_Tc ->{k_l} Sb_T_SOR + Sa_T_SOR-Tc_SOR_Tc	(* left-leak (rev) *))
		r_86 = k_l * X[13] * X[29] # (Sa_T_SOR + Sb_T_SOR-Tc_SOR_Tc ->{k_l} Sb_T_SOR + Sa_T_SOR-Tc_SOR_Tc	(* left-leak *))
		r_87 = k_l * X[24] * X[18] # (Sb_T_SOR + Sa_T_SOR-Tc_SOR_Tc ->{k_l} Sa_T_SOR + Sb_T_SOR-Tc_SOR_Tc	(* left-leak (rev) *))
		r_88 = k_l * X[13] * X[18] # (Sa_T_SOR + Sa_T_SOR-Tc_SOR_Tc ->{k_l} Sa_T_SOR + Sa_T_SOR-Tc_SOR_Tc	(* left-leak *))
		r_89 = k_l * X[13] * X[18] # (Sa_T_SOR + Sa_T_SOR-Tc_SOR_Tc ->{k_l} Sa_T_SOR + Sa_T_SOR-Tc_SOR_Tc	(* left-leak (rev) *))
		r_90 = k_l * X[39] * X[44] # (Ss10_T_Sf + Ss10_T_Sf-Tc_Ss10_Tc ->{k_l} Ss10_T_Sf + Ss10_T_Sf-Tc_Ss10_Tc	(* right-leak *))
		r_91 = k_l * X[39] * X[44] # (Ss10_T_Sf + Ss10_T_Sf-Tc_Ss10_Tc ->{k_l} Ss10_T_Sf + Ss10_T_Sf-Tc_Ss10_Tc	(* right-leak (rev) *))
		r_92 = k_l * X[39] * X[51] # (Ss10_T_Sf + Ss10_T_Sout-Tc_Ss10_Tc ->{k_l} Ss10_T_Sout + Ss10_T_Sf-Tc_Ss10_Tc	(* right-leak *))
		r_93 = k_l * X[45] * X[44] # (Ss10_T_Sout + Ss10_T_Sf-Tc_Ss10_Tc ->{k_l} Ss10_T_Sf + Ss10_T_Sout-Tc_Ss10_Tc	(* right-leak (rev) *))
		r_94 = k_l * X[45] * X[44] # (Ss10_T_Sout + Ss10_T_Sf-Tc_Ss10_Tc ->{k_l} Ss10_T_Sf + Ss10_T_Sout-Tc_Ss10_Tc	(* right-leak *))
		r_95 = k_l * X[39] * X[51] # (Ss10_T_Sf + Ss10_T_Sout-Tc_Ss10_Tc ->{k_l} Ss10_T_Sout + Ss10_T_Sf-Tc_Ss10_Tc	(* right-leak (rev) *))
		r_96 = k_l * X[45] * X[51] # (Ss10_T_Sout + Ss10_T_Sout-Tc_Ss10_Tc ->{k_l} Ss10_T_Sout + Ss10_T_Sout-Tc_Ss10_Tc	(* right-leak *))
		r_97 = k_l * X[45] * X[51] # (Ss10_T_Sout + Ss10_T_Sout-Tc_Ss10_Tc ->{k_l} Ss10_T_Sout + Ss10_T_Sout-Tc_Ss10_Tc	(* right-leak (rev) *))
		r_98 = k_l * X[0] * X[6] # (SOR_T_Ss10 + SOR_T_Ss10-Tc_Ss10_Tc ->{k_l} SOR_T_Ss10 + SOR_T_Ss10-Tc_Ss10_Tc	(* left-leak *))
		r_99 = k_l * X[0] * X[6] # (SOR_T_Ss10 + SOR_T_Ss10-Tc_Ss10_Tc ->{k_l} SOR_T_Ss10 + SOR_T_Ss10-Tc_Ss10_Tc	(* left-leak (rev) *))
		
		# species derivatives:
		dsp_0 = +r_0 -r_1 +r_2 -r_3 -r_4 +r_5 -r_6 +r_7 -r_8 -r_66 +r_67 -r_68 +r_69 -r_70 +r_71 -r_72 +r_73 -r_74 +r_75 -r_76 +r_77 -r_78 +r_79 -r_80 +r_80 -r_81 +r_81 -r_98 +r_98 -r_99 +r_99
		dsp_1 = +r_66 -r_67
		dsp_2 = +r_78 -r_79
		dsp_3 = +r_72 -r_73
		dsp_4 = +r_74 -r_75
		dsp_5 = -r_0 +r_1 -r_2 +r_3 -r_10 +r_11 -r_24 +r_25 -r_38 +r_39 -r_52 +r_53 -r_66 +r_67 -r_80 +r_80 -r_81 +r_81
		dsp_6 = +r_4 -r_5 +r_6 -r_7 -r_20 +r_21 -r_34 +r_35 -r_48 +r_49 -r_62 +r_63 -r_76 +r_77 -r_98 +r_98 -r_99 +r_99
		dsp_7 = +r_76 -r_77
		dsp_8 = +r_48 -r_49
		dsp_9 = +r_34 -r_35
		dsp_10 = +r_20 -r_21
		dsp_11 = +r_62 -r_63
		dsp_12 = +r_8
		dsp_13 = -r_2 +r_3 -r_38 +r_39 -r_40 +r_41 -r_42 +r_43 -r_44 +r_45 -r_46 +r_47 -r_48 +r_49 -r_50 +r_51 +r_84 -r_85 -r_86 +r_87 -r_88 +r_88 -r_89 +r_89
		dsp_14 = +r_38 -r_39
		dsp_15 = +r_50 -r_51
		dsp_16 = +r_44 -r_45
		dsp_17 = +r_46 -r_47
		dsp_18 = +r_2 -r_3 -r_14 +r_15 -r_28 +r_29 -r_42 +r_43 -r_56 +r_57 -r_70 +r_71 -r_84 +r_85 +r_86 -r_87 -r_88 +r_88 -r_89 +r_89
		dsp_19 = +r_70 -r_71
		dsp_20 = +r_42 -r_43
		dsp_21 = +r_28 -r_29
		dsp_22 = +r_14 -r_15
		dsp_23 = +r_56 -r_57
		dsp_24 = -r_0 +r_1 -r_24 +r_25 -r_26 +r_27 -r_28 +r_29 -r_30 +r_31 -r_32 +r_33 -r_34 +r_35 -r_36 +r_37 -r_82 +r_82 -r_83 +r_83 -r_84 +r_85 +r_86 -r_87
		dsp_25 = +r_24 -r_25
		dsp_26 = +r_36 -r_37
		dsp_27 = +r_30 -r_31
		dsp_28 = +r_32 -r_33
		dsp_29 = +r_0 -r_1 -r_12 +r_13 -r_26 +r_27 -r_40 +r_41 -r_54 +r_55 -r_68 +r_69 -r_82 +r_82 -r_83 +r_83 +r_84 -r_85 -r_86 +r_87
		dsp_30 = +r_68 -r_69
		dsp_31 = +r_40 -r_41
		dsp_32 = +r_26 -r_27
		dsp_33 = +r_12 -r_13
		dsp_34 = +r_54 -r_55
		dsp_35 = +r_9
		dsp_36 = -r_9
		dsp_37 = +r_8
		dsp_38 = -r_8 -r_22 +r_23 -r_36 +r_37 -r_50 +r_51 -r_64 +r_65 -r_78 +r_79
		dsp_39 = +r_4 -r_5 -r_10 +r_11 -r_12 +r_13 -r_14 +r_15 -r_16 +r_17 -r_18 +r_19 -r_20 +r_21 -r_22 +r_23 -r_90 +r_90 -r_91 +r_91 -r_92 +r_93 +r_94 -r_95
		dsp_40 = +r_10 -r_11
		dsp_41 = +r_22 -r_23
		dsp_42 = +r_16 -r_17
		dsp_43 = +r_18 -r_19
		dsp_44 = -r_4 +r_5 -r_16 +r_17 -r_30 +r_31 -r_44 +r_45 -r_58 +r_59 -r_72 +r_73 -r_90 +r_90 -r_91 +r_91 +r_92 -r_93 -r_94 +r_95
		dsp_45 = +r_6 -r_7 -r_9 -r_52 +r_53 -r_54 +r_55 -r_56 +r_57 -r_58 +r_59 -r_60 +r_61 -r_62 +r_63 -r_64 +r_65 +r_92 -r_93 -r_94 +r_95 -r_96 +r_96 -r_97 +r_97
		dsp_46 = +r_52 -r_53
		dsp_47 = +r_64 -r_65
		dsp_48 = +r_58 -r_59
		dsp_49 = +r_60 -r_61
		dsp_50 = +r_9
		dsp_51 = -r_6 +r_7 -r_18 +r_19 -r_32 +r_33 -r_46 +r_47 -r_60 +r_61 -r_74 +r_75 -r_92 +r_93 +r_94 -r_95 -r_96 +r_96 -r_97 +r_97
		
		return [ dsp_0, dsp_1, dsp_2, dsp_3, dsp_4, dsp_5, dsp_6, dsp_7, dsp_8, dsp_9, dsp_10, dsp_11, dsp_12, dsp_13, dsp_14, dsp_15, dsp_16, dsp_17, dsp_18, dsp_19, dsp_20, dsp_21, dsp_22, dsp_23, dsp_24, dsp_25, dsp_26, dsp_27, dsp_28, dsp_29, dsp_30, dsp_31, dsp_32, dsp_33, dsp_34, dsp_35, dsp_36, dsp_37, dsp_38, dsp_39, dsp_40, dsp_41, dsp_42, dsp_43, dsp_44, dsp_45, dsp_46, dsp_47, dsp_48, dsp_49, dsp_50, dsp_51 ]
	
	t = linspace(0, max_t, num = num_samples)
	values = integrate.odeint(dXdY, initial, t)
	
	ret = []
	
	# get the final states of the species asked for
	for specimen in get_final_states:
		ret.append(values[-1:,all_species.index(specimen)][0])
	return ret # return final states

#print analyze_system(12 * 60 * 60, 1000, [ SoutQ ])