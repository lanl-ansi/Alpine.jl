"""
// @brief  Script to generate the nonconvex instances with only bilinear terms
//		   Input arguments: numVariables, numInstances
// Reference: "Learning to Accelerate the Global Optimization of Quadratically-Constrained 
			   Quadratic Programs" by R. Kannan, H. Nagarajan, and D. Deka 
			   (Link: https://arxiv.org/abs/2301.00306)
// Contributor: @rohitkannan 
"""

import random
import json
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--numVariables', type=int, required=True)
parser.add_argument('--numInstances', type=int, required=True)
args = parser.parse_args()

numVariables = int(args.numVariables)
numInstances = int(args.numInstances)

instances_path = "./N_" + str(numVariables) + "/"
os.makedirs(instances_path)

# random seed for instances
if numVariables == 10:
	random.seed(4898)
elif numVariables == 20:
	random.seed(6736) 
elif numVariables == 50:
	random.seed(8197) 

np.random.seed(round(random.random()*10000)) 


start = 1  # starting instance number

numQuadraticConstraints = numVariables  # number of constraints with bilinear and linear terms
numLinearConstraints = round(numVariables/5)  # number of linear equality constraints

constraint_sparsity = 1.0  # (actually opposite of sparsity) fraction of possible bilinear and linear terms occurring in each constraint
maxNumBilinearTerms = round(numVariables*(numVariables-1)/2)
numBilinearTerms = min(5*numVariables,maxNumBilinearTerms)  # number of bilinear terms in the formulation

frac_quad_con_pert = 0.2  # fraction of quadratic constraints perturbed
numPertFactors = 3  # number of perturbation factors per perturbed quadratic constraint
pert_frac = 0.5  # perturbation fraction around nominal instance


# to sort list of lists
def Sort(sub_li):
    return(sorted(sub_li, key = lambda x: (x[0],x[1])))    


# GENERATE COMMON BILINEAR TERMS FOR ALL INSTANCES
possible_bilinear_terms = [[i,j] for i in range(numVariables) for j in range(i+1,numVariables)]
random.shuffle(possible_bilinear_terms)
bilinear_terms = Sort(possible_bilinear_terms[:numBilinearTerms])


# GENERATE BASE COEFFICIENTS
q_bilinear = []  # bilinear coefficients
q_bilinear_ind = []  # bilinear coefficient indices

# OBJECTIVE COEFFICIENTS
q_bilinear.append([2*random.random()-1 for j in range(numBilinearTerms)])  # all bilinear terms participate in the objective (sparsity not applied here)
q_bilinear_ind.append(list(np.arange(numBilinearTerms)))

# CONSTRAINT COEFFICIENTS
numQuadNonZero = round(numBilinearTerms*constraint_sparsity)
for i in range(numQuadraticConstraints):
	q_bilinear.append([2*random.random()-1 for j in range(numQuadNonZero)])
	tmp_ind = list(np.random.choice(numBilinearTerms, numQuadNonZero, replace=False))
	tmp_ind.sort()
	q_bilinear_ind.append(tmp_ind)


q_linear = []  # linear coefficients
q_linear_ind = []  # linear coefficient indices

# OBJECTIVE COEFFICIENTS
q_linear.append([2*random.random()-1 for j in range(numVariables)])  # all linear terms participate in the objective (sparsity not applied here)
q_linear_ind.append(list(np.arange(numVariables)))

# CONSTRAINT COEFFICIENTS
numLinNonZero = round(numVariables*constraint_sparsity)
for i in range(numQuadraticConstraints+numLinearConstraints):
	if i < numQuadraticConstraints:  # only a fraction of linear terms appear in quadratic constraints
		q_linear.append([2*random.random()-1 for j in range(numLinNonZero)])
		tmp_ind = list(np.random.choice(numVariables, numLinNonZero, replace=False))
	else:  # all linear terms appear in linear equality constraints
		q_linear.append([2*random.random()-1 for j in range(numVariables)])
		tmp_ind = list(np.arange(numVariables))
	tmp_ind.sort()
	q_linear_ind.append(tmp_ind)


# GENERATE BASE RHS COEFFICIENTS
quadratic_scaling = 100
rhs_quadratic = [quadratic_scaling*random.random() for i in range(numQuadraticConstraints)]

linear_scaling = 1
rhs_linear = [linear_scaling*(2*random.random()-1) for i in range(numLinearConstraints)]


# RE-SCALE PROBLEM SO THAT CONSTRAINTS RHS ARE EQUAL TO ONE
for i in range(numQuadraticConstraints):
	for j in range(numQuadNonZero):
		q_bilinear[i+1][j] /= rhs_quadratic[i]
	for j in range(numLinNonZero):
		q_linear[i+1][j] /= rhs_quadratic[i]


for i in range(numLinearConstraints):
	for j in range(numVariables):
		q_linear[i+1+numQuadraticConstraints][j] /= rhs_linear[i]


# CONSTRUCT PERTURBATION FACTORS FOR THE OBJECTIVE AND QUADRATIC CONSTRAINT COEFFICIENTS
q_bilinear_obj_factors = np.zeros((numPertFactors,numBilinearTerms))
q_linear_obj_factors = np.zeros((numPertFactors,numVariables))
for j in range(numBilinearTerms):
	for k in range(numPertFactors):
		q_bilinear_obj_factors[k,j] = pert_frac*random.random()*q_bilinear[0][j]
for j in range(numVariables):
	for k in range(numPertFactors):
		q_linear_obj_factors[k,j] = pert_frac*random.random()*q_linear[0][j]

numQuadConPert = round(numQuadraticConstraints*frac_quad_con_pert)
q_bilinear_con_factors = np.zeros((numPertFactors,numQuadConPert,numQuadNonZero))
q_linear_con_factors = np.zeros((numPertFactors,numQuadConPert,numLinNonZero))
for i in range(numQuadConPert):
	for j in range(numQuadNonZero):
		for k in range(numPertFactors):
			q_bilinear_con_factors[k,i,j] = pert_frac*random.random()*q_bilinear[i+1][j]
	for j in range(numLinNonZero):
		for k in range(numPertFactors):
			q_linear_con_factors[k,i,j] = pert_frac*random.random()*q_linear[i+1][j]



# GENERATE PERTURBED INSTANCES
for inst in range(start,numInstances+start):
	output_file = instances_path + "bilinear_v" + str(numVariables) + "_b" + str(numBilinearTerms) + "_s" + str(round(constraint_sparsity*100)) + "_" + str(inst) + ".jl"
	json_file = instances_path + "bilinear_v" + str(numVariables) + "_b" + str(numBilinearTerms) + "_s" + str(round(constraint_sparsity*100)) + "_" + str(inst) + ".json"

	json_dict = {}
	json_dict["numVariables"] = numVariables
	json_dict["numQuadraticConstraints"] = numQuadraticConstraints
	json_dict["numLinearConstraints"] = numLinearConstraints
	json_dict["constraint_sparsity"] = constraint_sparsity

	json_dict["numBilinearTerms"] = numBilinearTerms
	json_dict["bilinear_terms"] = bilinear_terms

	json_dict["pert_frac"] = pert_frac
	json_dict["numPertFactors"] = numPertFactors
	json_dict["numQuadConPert"] = numQuadConPert
	json_dict["bilinear_obj_factors"] = q_bilinear_obj_factors.tolist()
	json_dict["linear_obj_factors"] = q_linear_obj_factors.tolist()
	json_dict["bilinear_con_factors"] = q_bilinear_con_factors.tolist()
	json_dict["linear_con_factors"] = q_linear_con_factors.tolist()


# GENERATE PERTURBATION FACTORS
	theta = 2*(np.random.rand(numPertFactors,numQuadConPert+1) - 0.5)

	json_dict["theta"] = theta.tolist()

# STORE SPARSITY PATTERN AND INITIALIZE DATA
	f = open(output_file, "w")

	f.write("# ----- Bilinear Terms ----- #\n")
	f.write("numBilinearTerms = " + str(numBilinearTerms) + "\n")
	f.write("bilinear_terms = [")
	count_bilinear = 0
	for term in bilinear_terms:
		f.write("[" + str(term[0]+1) + "," + str(term[1]+1) + "]")
		count_bilinear += 1
		if count_bilinear < numBilinearTerms:
			f.write(", ")
		else:
			f.write("] \n\n")
	
	f.write("numQuadNonZero = " + str(numQuadNonZero) + "\n")
	f.write("numLinNonZero = " + str(numLinNonZero) + "\n\n")

	f.write("# ----- Data Matrices ----- #\n")
	f.write("bilinear_obj_coeffs = zeros(Float64," + str(numBilinearTerms) + ") \n")
	f.write("bilinear_con_coeffs = zeros(Float64," + str(numQuadraticConstraints) + "," + str(numQuadNonZero) + ") \n")
	f.write("bilinear_con_coeff_indices = zeros(Int64," + str(numQuadraticConstraints) + "," + str(numQuadNonZero) + ") \n")
	f.write("linear_obj_coeffs = zeros(Float64," + str(numVariables) + ") \n")
	f.write("linear_con_coeffs_1 = zeros(Float64," + str(numQuadraticConstraints) + "," + str(numLinNonZero) + ") \n")
	f.write("linear_con_coeff_indices_1 = zeros(Int64," + str(numQuadraticConstraints) + "," + str(numLinNonZero) + ") \n")
	f.write("linear_con_coeffs_2 = zeros(Float64," + str(numLinearConstraints) + "," + str(numVariables) + ") \n")
	f.write("linear_con_coeff_indices_2 = zeros(Int64," + str(numLinearConstraints) + "," + str(numVariables) + ") \n\n")

	f.write("# ----- Variables ----- #\n")
	f.write("@variable(m, 0 <= x[1:" + str(numVariables) + "] <= 1)\n")
	f.write("@variable(m, 0 <= w[1:" + str(numBilinearTerms) + "] <= 1)\n\n")

	f.write("# ----- Bilinear Equations ----- #\n")
	f.write("@NLconstraint(m, [i=1:" + str(numBilinearTerms) + "], w[i] == x[bilinear_terms[i][1]]*x[bilinear_terms[i][2]])\n\n\n")

	json_dict["x_lower"] = [0 for i in range(numVariables)]
	json_dict["x_upper"] = [1 for i in range(numVariables)]


# OBJECTIVE
	# tmp variables for storing generated data
	Q = [0.]*numBilinearTerms
	q = [0.]*numVariables

	count_quad_con = 0
	count_lin_con = 0
	for i in range(numVariables):
		q[i] = q_linear[0][i]
		for k in range(numPertFactors):
			q[i] += q_linear_obj_factors[k,i]*theta[k,0]
	for i in range(numBilinearTerms):
		Q[i] = q_bilinear[0][i]
		for k in range(numPertFactors):
			Q[i] += q_bilinear_obj_factors[k,i]*theta[k,0]

	json_dict["obj_bilinear_coeffs"] = Q
	json_dict["obj_linear_coeffs"] = q

	f.write("# ----- Objective ----- #\n")
	f.write("bilinear_obj_coeffs = " + str(Q) + "\n")
	f.write("linear_obj_coeffs = " + str(q) + "\n")

	f.write("\n")
	f.write("@objective(m, Min, sum(bilinear_obj_coeffs[i]*w[i] for i = 1:" + str(numBilinearTerms) + ") + sum(linear_obj_coeffs[i]*x[i] for i = 1:" + str(numVariables) + ")) \n\n\n")


	bilinear_obj_matrix = np.zeros((numVariables, numVariables))
	for i in range(numBilinearTerms):
		bilinear_obj_matrix[bilinear_terms[i][0]][bilinear_terms[i][1]] = Q[i]/2.0
		bilinear_obj_matrix[bilinear_terms[i][1]][bilinear_terms[i][0]] = Q[i]/2.0

	json_dict["obj_bilinear_eigvals"] = np.linalg.eigvals(bilinear_obj_matrix).tolist()


# QUADRATIC CONSTRAINTS
	f.write("# ----- Constraints ----- #\n")

	json_dict["quadratic_con_rhs"] = [1.0 for i in range(numQuadraticConstraints)]

	bilinear_con_coeffs = []
	bilinear_con_coeff_indices = []
	linear_con_coeffs = []
	linear_con_coeff_indices = []
	con_bilinear_eigvals = np.empty([numQuadraticConstraints, numVariables], dtype=float)

	for k in range(numQuadraticConstraints):
		Q = [0.]*numQuadNonZero
		Q_ind = [i+1 for i in q_bilinear_ind[count_quad_con+1]]
		q = [0.]*numLinNonZero
		q_ind = [i+1 for i in q_linear_ind[count_lin_con+1]]

		count_quad_con += 1
		count_lin_con += 1
		for i in range(numLinNonZero):
			q[i] = q_linear[count_lin_con][i]
			if k < numQuadConPert:
				for j in range(numPertFactors):
					q[i] += q_linear_con_factors[j,count_lin_con-1,i]*theta[j,count_lin_con]
		for i in range(numQuadNonZero):
			Q[i] = q_bilinear[count_quad_con][i]
			if k < numQuadConPert:
				for j in range(numPertFactors):
					Q[i] += q_bilinear_con_factors[j,count_quad_con-1,i]*theta[j,count_quad_con]

		bilinear_con_coeffs.append(Q)
		bilinear_con_coeff_indices.append([float(i) for i in q_bilinear_ind[count_quad_con]])
		linear_con_coeffs.append(q)
		linear_con_coeff_indices.append([float(i) for i in q_linear_ind[count_lin_con]])

		f.write("bilinear_con_coeffs[" + str(count_quad_con) + ",:] = " + str(Q) + "\n")
		f.write("bilinear_con_coeff_indices[" + str(count_quad_con) + ",:] = " + str(Q_ind) + "\n")
		f.write("linear_con_coeffs_1[" + str(count_lin_con) + ",:] = " + str(q) + "\n")
		f.write("linear_con_coeff_indices_1[" + str(count_lin_con) + ",:] = " + str(q_ind) + "\n")
		f.write("\n")

		f.write("@constraint(m, sum(bilinear_con_coeffs[" + str(count_quad_con) + ",i]*w[bilinear_con_coeff_indices[" + str(count_quad_con) + ",i]] for i = 1:" + str(numQuadNonZero) + ") + sum(linear_con_coeffs_1[" + str(count_lin_con) + ",i]*x[linear_con_coeff_indices_1[" + str(count_lin_con) + ",i]] for i = 1:" + str(numLinNonZero) + ") <= 1) \n\n\n")


		bilinear_con_matrix = np.zeros((numVariables, numVariables))
		for i in range(numQuadNonZero):
			ind = q_bilinear_ind[count_quad_con][i]
			ind1 = bilinear_terms[ind][0]
			ind2 = bilinear_terms[ind][1]
			bilinear_con_matrix[ind1][ind2] = Q[i]/2.0
			bilinear_con_matrix[ind2][ind1] = Q[i]/2.0

		con_bilinear_eigvals[k] = np.linalg.eigvals(bilinear_con_matrix)

	
	json_dict["con_bilinear_eigvals"] = con_bilinear_eigvals.tolist()


# LINEAR CONSTRAINTS
	json_dict["linear_con_rhs"] = [1.0 for i in range(numLinearConstraints)]

	for k in range(numLinearConstraints):
		q = [0.]*numVariables
		q_ind = [i+1 for i in q_linear_ind[count_lin_con+1]]

		count_lin_con += 1
		for i in range(numVariables):
			q[i] = q_linear[count_lin_con][i]

		linear_con_coeffs.append(q)
		linear_con_coeff_indices.append([float(i) for i in q_linear_ind[count_lin_con]])

		f.write("linear_con_coeffs_2[" + str(count_lin_con - numQuadraticConstraints) + ",:] = " + str(q) + "\n")
		f.write("linear_con_coeff_indices_2[" + str(count_lin_con - numQuadraticConstraints) + ",:] = " + str(q_ind) + "\n")
		f.write("\n")

		f.write("@constraint(m, sum(linear_con_coeffs_2[" + str(count_lin_con - numQuadraticConstraints) + ",i]*x[linear_con_coeff_indices_2[" + str(count_lin_con - numQuadraticConstraints) + ",i]] for i = 1:" + str(numVariables) + ") == 1) \n\n\n")
		
	f.close()


	json_dict["con_bilinear_coeffs"] = bilinear_con_coeffs
	json_dict["con_bilinear_coeff_indices"] = bilinear_con_coeff_indices
	json_dict["con_linear_coeffs"] = linear_con_coeffs
	json_dict["con_linear_coeff_indices"] = linear_con_coeff_indices

	with open(json_file, 'w', encoding='utf-8') as outF:
		json.dump(json_dict, outF, ensure_ascii=False, indent=4)
