from scipy import integrate
from numpy import linspace

def analyze_system(max_t, num_samples, A_amount, get_final_states):
	""" Returns the final state of interesting species for the system "NOT gate" """
	all_species = ['SNOT_Ss10-Tc_SNOTc_Ss10c_Tc', 'SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11', 'Sa_T_SNOT', 'Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc', 'Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11', 'Sa_T_SNOT-Ss11_T_Sf-Tc_Ss11_Tc', 'Sa_T_SNOT-Ss11_T_Sout-Tc_Ss11_Tc', 'SoutQ', 'SoutQ-Tc_SoutcF', 'Ss10_T_Ss11', 'Ss10_T_Ss11-Ss11_T_Sf-Tc_Ss11_Tc', 'Ss10_T_Ss11-Ss11_T_Sout-Tc_Ss11_Tc', 'Ss10_T_Ss11-Tc_Ss11_Tc', 'Ss10_T_Ss11-Tc_Ss11_Tc-Sa_T_SNOT', 'Ss10_T_Ss11-Tc_Ss11_Tc-Ss10_T_Ss11', 'Ss10_T_Ss11-Tc_Ss11_Tc-Ss11_T_Sf', 'Ss10_T_Ss11-Tc_Ss11_Tc-Ss11_T_Sout', 'Ss11_T_Sf', 'Ss11_T_Sf-Ss11_T_Sf-Tc_Ss11_Tc', 'Ss11_T_Sf-Ss11_T_Sout-Tc_Ss11_Tc', 'Ss11_T_Sf-Tc_Ss11_Tc', 'Ss11_T_Sout', 'Ss11_T_Sout-Ss11_T_Sf-Tc_Ss11_Tc', 'Ss11_T_Sout-Ss11_T_Sout-Tc_Ss11_Tc', 'Ss11_T_Sout-Tc_SoutcF', 'Ss11_T_Sout-Tc_Ss11_Tc']
	num_species = len(all_species)
	
	initial = [ 0.0 ] * num_species # initial conditions
	initial[8] = 150 # SoutQ-Tc_SoutcF
	initial[2] = A_amount * 100 # Sa_T_SNOT
	initial[25] = 1.0 * 100 # Ss11_T_Sout-Tc_Ss11_Tc
	initial[9] = 0.5 * 100 # Ss10_T_Ss11
	initial[17] = 2.0 * 100 # Ss11_T_Sf
	initial[0] = 10.0 * 100 # SNOT_Ss10-Tc_SNOTc_Ss10c_Tc
	
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
		r_0 = k_s * X[9] * X[25] # (Ss10_T_Ss11 + Ss11_T_Sout-Tc_Ss11_Tc ->{k_s} Ss10_T_Ss11-Tc_Ss11_Tc + Ss11_T_Sout	(* see-saw *))
		r_1 = k_s * X[12] * X[21] # (Ss10_T_Ss11-Tc_Ss11_Tc + Ss11_T_Sout ->{k_s} Ss10_T_Ss11 + Ss11_T_Sout-Tc_Ss11_Tc	(* see-saw (rev) *))
		r_2 = k_s * X[9] * X[20] # (Ss10_T_Ss11 + Ss11_T_Sf-Tc_Ss11_Tc ->{k_s} Ss10_T_Ss11-Tc_Ss11_Tc + Ss11_T_Sf	(* see-saw *))
		r_3 = k_s * X[12] * X[17] # (Ss10_T_Ss11-Tc_Ss11_Tc + Ss11_T_Sf ->{k_s} Ss10_T_Ss11 + Ss11_T_Sf-Tc_Ss11_Tc	(* see-saw (rev) *))
		r_4 = k_s * X[2] * X[0] # (Sa_T_SNOT + SNOT_Ss10-Tc_SNOTc_Ss10c_Tc ->{k_s} Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc	(* see-saw *))
		r_5 = k_s * X[3] # (Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc ->{k_s} Sa_T_SNOT + SNOT_Ss10-Tc_SNOTc_Ss10c_Tc	(* see-saw (rev) *))
		r_6 = k_s * X[3] * X[9] # (Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc + Ss10_T_Ss11 ->{k_s} Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11	(* see-saw (irreversible cooperative hybridization) *))
		r_7 = k_s * X[9] * X[0] # (Ss10_T_Ss11 + SNOT_Ss10-Tc_SNOTc_Ss10c_Tc ->{k_s} SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11	(* see-saw *))
		r_8 = k_s * X[1] # (SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11 ->{k_s} Ss10_T_Ss11 + SNOT_Ss10-Tc_SNOTc_Ss10c_Tc	(* see-saw (rev) *))
		r_9 = k_s * X[1] * X[2] # (SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11 + Sa_T_SNOT ->{k_s} Sa_T_SNOT-SNOT_Ss10-Tc_SNOTc_Ss10c_Tc-Ss10_T_Ss11	(* see-saw (irreversible cooperative hybridization) *))
		r_10 = k_s * X[21] * X[8] # (Ss11_T_Sout + SoutQ-Tc_SoutcF ->{k_s} Ss11_T_Sout-Tc_SoutcF + SoutQ	(* reporting *))
		r_11 = k_f * X[9] * X[25] # (Ss10_T_Ss11 + Ss11_T_Sout-Tc_Ss11_Tc ->{k_f} Ss10_T_Ss11-Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_12 = k_rf * X[11] # (Ss10_T_Ss11-Ss11_T_Sout-Tc_Ss11_Tc ->{k_rf} Ss10_T_Ss11 + Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_13 = k_f * X[9] * X[20] # (Ss10_T_Ss11 + Ss11_T_Sf-Tc_Ss11_Tc ->{k_f} Ss10_T_Ss11-Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_14 = k_rf * X[10] # (Ss10_T_Ss11-Ss11_T_Sf-Tc_Ss11_Tc ->{k_rf} Ss10_T_Ss11 + Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_15 = k_f * X[9] * X[12] # (Ss10_T_Ss11 + Ss10_T_Ss11-Tc_Ss11_Tc ->{k_f} Ss10_T_Ss11-Tc_Ss11_Tc-Ss10_T_Ss11	(* universal toehold leak: from right *))
		r_16 = k_rf * X[14] # (Ss10_T_Ss11-Tc_Ss11_Tc-Ss10_T_Ss11 ->{k_rf} Ss10_T_Ss11 + Ss10_T_Ss11-Tc_Ss11_Tc	(* universal toehold leak: from right (rev) *))
		r_17 = k_f * X[21] * X[25] # (Ss11_T_Sout + Ss11_T_Sout-Tc_Ss11_Tc ->{k_f} Ss11_T_Sout-Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_18 = k_rf * X[23] # (Ss11_T_Sout-Ss11_T_Sout-Tc_Ss11_Tc ->{k_rf} Ss11_T_Sout + Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_19 = k_f * X[21] * X[20] # (Ss11_T_Sout + Ss11_T_Sf-Tc_Ss11_Tc ->{k_f} Ss11_T_Sout-Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_20 = k_rf * X[22] # (Ss11_T_Sout-Ss11_T_Sf-Tc_Ss11_Tc ->{k_rf} Ss11_T_Sout + Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_21 = k_f * X[21] * X[12] # (Ss11_T_Sout + Ss10_T_Ss11-Tc_Ss11_Tc ->{k_f} Ss10_T_Ss11-Tc_Ss11_Tc-Ss11_T_Sout	(* universal toehold leak: from right *))
		r_22 = k_rf * X[16] # (Ss10_T_Ss11-Tc_Ss11_Tc-Ss11_T_Sout ->{k_rf} Ss11_T_Sout + Ss10_T_Ss11-Tc_Ss11_Tc	(* universal toehold leak: from right (rev) *))
		r_23 = k_f * X[17] * X[25] # (Ss11_T_Sf + Ss11_T_Sout-Tc_Ss11_Tc ->{k_f} Ss11_T_Sf-Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_24 = k_rf * X[19] # (Ss11_T_Sf-Ss11_T_Sout-Tc_Ss11_Tc ->{k_rf} Ss11_T_Sf + Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_25 = k_f * X[17] * X[20] # (Ss11_T_Sf + Ss11_T_Sf-Tc_Ss11_Tc ->{k_f} Ss11_T_Sf-Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_26 = k_rf * X[18] # (Ss11_T_Sf-Ss11_T_Sf-Tc_Ss11_Tc ->{k_rf} Ss11_T_Sf + Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_27 = k_f * X[17] * X[12] # (Ss11_T_Sf + Ss10_T_Ss11-Tc_Ss11_Tc ->{k_f} Ss10_T_Ss11-Tc_Ss11_Tc-Ss11_T_Sf	(* universal toehold leak: from right *))
		r_28 = k_rf * X[15] # (Ss10_T_Ss11-Tc_Ss11_Tc-Ss11_T_Sf ->{k_rf} Ss11_T_Sf + Ss10_T_Ss11-Tc_Ss11_Tc	(* universal toehold leak: from right (rev) *))
		r_29 = k_f * X[2] * X[25] # (Sa_T_SNOT + Ss11_T_Sout-Tc_Ss11_Tc ->{k_f} Sa_T_SNOT-Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_30 = k_rf * X[6] # (Sa_T_SNOT-Ss11_T_Sout-Tc_Ss11_Tc ->{k_rf} Sa_T_SNOT + Ss11_T_Sout-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_31 = k_f * X[2] * X[20] # (Sa_T_SNOT + Ss11_T_Sf-Tc_Ss11_Tc ->{k_f} Sa_T_SNOT-Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left *))
		r_32 = k_rf * X[5] # (Sa_T_SNOT-Ss11_T_Sf-Tc_Ss11_Tc ->{k_rf} Sa_T_SNOT + Ss11_T_Sf-Tc_Ss11_Tc	(* universal toehold leak: from left (rev) *))
		r_33 = k_f * X[2] * X[12] # (Sa_T_SNOT + Ss10_T_Ss11-Tc_Ss11_Tc ->{k_f} Ss10_T_Ss11-Tc_Ss11_Tc-Sa_T_SNOT	(* universal toehold leak: from right *))
		r_34 = k_rf * X[13] # (Ss10_T_Ss11-Tc_Ss11_Tc-Sa_T_SNOT ->{k_rf} Sa_T_SNOT + Ss10_T_Ss11-Tc_Ss11_Tc	(* universal toehold leak: from right (rev) *))
		r_35 = k_l * X[21] * X[25] # (Ss11_T_Sout + Ss11_T_Sout-Tc_Ss11_Tc ->{k_l} Ss11_T_Sout + Ss11_T_Sout-Tc_Ss11_Tc	(* right-leak *))
		r_36 = k_l * X[21] * X[25] # (Ss11_T_Sout + Ss11_T_Sout-Tc_Ss11_Tc ->{k_l} Ss11_T_Sout + Ss11_T_Sout-Tc_Ss11_Tc	(* right-leak (rev) *))
		r_37 = k_l * X[21] * X[20] # (Ss11_T_Sout + Ss11_T_Sf-Tc_Ss11_Tc ->{k_l} Ss11_T_Sf + Ss11_T_Sout-Tc_Ss11_Tc	(* right-leak *))
		r_38 = k_l * X[17] * X[25] # (Ss11_T_Sf + Ss11_T_Sout-Tc_Ss11_Tc ->{k_l} Ss11_T_Sout + Ss11_T_Sf-Tc_Ss11_Tc	(* right-leak (rev) *))
		r_39 = k_l * X[17] * X[25] # (Ss11_T_Sf + Ss11_T_Sout-Tc_Ss11_Tc ->{k_l} Ss11_T_Sout + Ss11_T_Sf-Tc_Ss11_Tc	(* right-leak *))
		r_40 = k_l * X[21] * X[20] # (Ss11_T_Sout + Ss11_T_Sf-Tc_Ss11_Tc ->{k_l} Ss11_T_Sf + Ss11_T_Sout-Tc_Ss11_Tc	(* right-leak (rev) *))
		r_41 = k_l * X[17] * X[20] # (Ss11_T_Sf + Ss11_T_Sf-Tc_Ss11_Tc ->{k_l} Ss11_T_Sf + Ss11_T_Sf-Tc_Ss11_Tc	(* right-leak *))
		r_42 = k_l * X[17] * X[20] # (Ss11_T_Sf + Ss11_T_Sf-Tc_Ss11_Tc ->{k_l} Ss11_T_Sf + Ss11_T_Sf-Tc_Ss11_Tc	(* right-leak (rev) *))
		r_43 = k_l * X[9] * X[12] # (Ss10_T_Ss11 + Ss10_T_Ss11-Tc_Ss11_Tc ->{k_l} Ss10_T_Ss11 + Ss10_T_Ss11-Tc_Ss11_Tc	(* left-leak *))
		r_44 = k_l * X[9] * X[12] # (Ss10_T_Ss11 + Ss10_T_Ss11-Tc_Ss11_Tc ->{k_l} Ss10_T_Ss11 + Ss10_T_Ss11-Tc_Ss11_Tc	(* left-leak (rev) *))
		
		# species derivatives:
		dsp_0 = -r_4 +r_5 -r_7 +r_8
		dsp_1 = +r_7 -r_8 -r_9
		dsp_2 = -r_4 +r_5 -r_9 -r_29 +r_30 -r_31 +r_32 -r_33 +r_34
		dsp_3 = +r_4 -r_5 -r_6
		dsp_4 = +r_6 +r_9
		dsp_5 = +r_31 -r_32
		dsp_6 = +r_29 -r_30
		dsp_7 = +r_10
		dsp_8 = -r_10
		dsp_9 = -r_0 +r_1 -r_2 +r_3 -r_6 -r_7 +r_8 -r_11 +r_12 -r_13 +r_14 -r_15 +r_16 -r_43 +r_43 -r_44 +r_44
		dsp_10 = +r_13 -r_14
		dsp_11 = +r_11 -r_12
		dsp_12 = +r_0 -r_1 +r_2 -r_3 -r_15 +r_16 -r_21 +r_22 -r_27 +r_28 -r_33 +r_34 -r_43 +r_43 -r_44 +r_44
		dsp_13 = +r_33 -r_34
		dsp_14 = +r_15 -r_16
		dsp_15 = +r_27 -r_28
		dsp_16 = +r_21 -r_22
		dsp_17 = +r_2 -r_3 -r_23 +r_24 -r_25 +r_26 -r_27 +r_28 +r_37 -r_38 -r_39 +r_40 -r_41 +r_41 -r_42 +r_42
		dsp_18 = +r_25 -r_26
		dsp_19 = +r_23 -r_24
		dsp_20 = -r_2 +r_3 -r_13 +r_14 -r_19 +r_20 -r_25 +r_26 -r_31 +r_32 -r_37 +r_38 +r_39 -r_40 -r_41 +r_41 -r_42 +r_42
		dsp_21 = +r_0 -r_1 -r_10 -r_17 +r_18 -r_19 +r_20 -r_21 +r_22 -r_35 +r_35 -r_36 +r_36 -r_37 +r_38 +r_39 -r_40
		dsp_22 = +r_19 -r_20
		dsp_23 = +r_17 -r_18
		dsp_24 = +r_10
		dsp_25 = -r_0 +r_1 -r_11 +r_12 -r_17 +r_18 -r_23 +r_24 -r_29 +r_30 -r_35 +r_35 -r_36 +r_36 +r_37 -r_38 -r_39 +r_40
		
		return [ dsp_0, dsp_1, dsp_2, dsp_3, dsp_4, dsp_5, dsp_6, dsp_7, dsp_8, dsp_9, dsp_10, dsp_11, dsp_12, dsp_13, dsp_14, dsp_15, dsp_16, dsp_17, dsp_18, dsp_19, dsp_20, dsp_21, dsp_22, dsp_23, dsp_24, dsp_25 ]
	
	t = linspace(0, max_t, num = num_samples)
	values = integrate.odeint(dXdY, initial, t)
	
	ret = []
	
	# get the final states of the species asked for
	for specimen in get_final_states:
		ret.append(values[-1:,all_species.index(specimen)][0])
	return ret # return final states

#print analyze_system(12 * 60 * 60, 1000, [ SoutQ ])
