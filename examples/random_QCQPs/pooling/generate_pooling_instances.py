"""
// @brief  Script to generate the nonconvex pooling instances based on Luedtke et al. (SIOPT vol. 30 (2020), pp. 1582–1609)
//		   Input arguments: numInstances
// Reference: "Learning to Accelerate the Global Optimization of Quadratically-Constrained 
			   Quadratic Programs" by R. Kannan, H. Nagarajan, and D. Deka 
			   (Link: https://arxiv.org/abs/2301.00306)
// Contributor: @rohitkannan 
"""

import numpy as np
import random
import json
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--numInstances', type=int, required=True)
args = parser.parse_args()

numInstances = int(args.numInstances)

instances_path = "./pooling_instances/"
os.mkdir(instances_path)

# random seed for instances
random.seed(8446) 
np.random.seed(round(random.random()*10000)) 


def flatten(t):
	return [item for sublist in t for item in sublist]


# parameters defining the instance
num_copies = 15
num_added_edges = 150
num_specs = 1


num_orig_inputs = 3
num_orig_pools = 1
num_orig_outputs = 2


# fraction by which to perturb the nominal input qualities
pert_frac = 0.2


# scaling factors for parameters
capacities_scaling = 100.0
costs_scaling = 2.0



input_file = "haverly_" + str(num_copies) + "_addedges_" + str(num_added_edges) + "_attr_" + str(num_specs-1) + "_1.dat"



# first map the vertex labels to integers
inputs_dict = dict()
outputs_dict = dict()
pools_dict = dict()

for c in range(num_copies):
	for i in range(num_orig_inputs):
		if c < 10 and num_copies > 10:
			tmp = "h0" + str(c) + "_i" + str(i+1)
		else:
			tmp = "h" + str(c) + "_i" + str(i+1)
		inputs_dict[tmp] = c*num_orig_inputs + i

	for j in range(num_orig_outputs):
		if c < 10 and num_copies > 10:
			tmp = "h0" + str(c) + "_j" + str(j+1)
		else:
			tmp = "h" + str(c) + "_j" + str(j+1)
		outputs_dict[tmp] = c*num_orig_outputs + j

	for l in range(num_orig_pools):
		if c < 10 and num_copies > 10:
			tmp = "h0" + str(c) + "_l" + str(l+1)
		else:
			tmp = "h" + str(c) + "_l" + str(l+1)
		pools_dict[tmp] = c*num_orig_pools + l


# compute actual number of inputs, pools, and outputs
num_inputs = num_orig_inputs*num_copies
num_pools = num_orig_pools*num_copies
num_outputs = num_orig_outputs*num_copies


# extract the relevant lines for each input parameter
file1 = open(input_file, 'r')
Lines = file1.readlines()

start_edges = False
start_costs = False
start_capacities = False
start_specs = False
start_limits = False

lines_edges = []
lines_costs = []
lines_capacities = []
lines_specs = []
lines_limits = []

for line in Lines:
	if start_edges and line.startswith('/;'):
		start_edges = False
	if start_costs and line.startswith('/;'):
		start_costs = False
	if start_capacities and line.startswith('/;'):
		start_capacities = False
	if start_specs and line.startswith('/;'):
		start_specs = False
	if start_limits and line.startswith('/;'):
		start_limits = False

	if start_edges:
		lines_edges.append(line)
	if start_costs:
		lines_costs.append(line)
	if start_capacities:
		lines_capacities.append(line)
	if start_specs:
		lines_specs.append(line)
	if start_limits:
		lines_limits.append(line)

	if line.startswith('set A /'):
		start_edges = True
	if line.startswith('parameter cost(V,V) /'):
		start_costs = True
	if line.startswith('parameter C(V) /'):
		start_capacities = True
	if line.startswith('parameter lambda(K,V) /'):
		start_specs = True
	if line.startswith('parameter overbeta(K,V) /'):
		start_limits = True


# first, extract capacities
cap_i = [0.0 for i in range(num_inputs)]
cap_j = [0.0 for j in range(num_outputs)]
cap_l = [0.0 for l in range(num_pools)]

for line in lines_capacities:
	if num_copies > 10:
		id = line[4:10]
		val = float(line[11:])
	else:
		id = line[4:9]
		val = float(line[10:])

	if id in inputs_dict:
		cap_i[inputs_dict[id]] = val/capacities_scaling
	if id in outputs_dict:
		cap_j[outputs_dict[id]] = val/capacities_scaling
	if id in pools_dict:
		cap_l[pools_dict[id]] = val/capacities_scaling


# next, extract costs
cost_ij_rows = []
cost_ij_cols = []
cost_ij_vals = []

cost_il_rows = []
cost_il_cols = []
cost_il_vals = []

cost_lj_rows = []
cost_lj_cols = []
cost_lj_vals = []

for line in lines_costs:
	ids = line.rstrip().split('.')
	if num_copies > 10:
		id1 = ids[0][-6:]
		id2 = ids[1][:6]
		val = float(line[18:])
	else:
		id1 = ids[0][-5:]
		id2 = ids[1][:5]
		val = float(line[16:])

	if id1 in inputs_dict and id2 in outputs_dict:
		cost_ij_rows.append(inputs_dict[id1])
		cost_ij_cols.append(outputs_dict[id2])
		cost_ij_vals.append(val/costs_scaling)

	if id1 in inputs_dict and id2 in pools_dict:
		cost_il_rows.append(inputs_dict[id1])
		cost_il_cols.append(pools_dict[id2])
		cost_il_vals.append(val/costs_scaling)

	if id1 in pools_dict and id2 in outputs_dict:
		cost_lj_rows.append(pools_dict[id1])
		cost_lj_cols.append(outputs_dict[id2])
		cost_lj_vals.append(val/costs_scaling)


# next, extract the edges
arcs_ij = []
arcs_il = []
arcs_lj = []

for line in lines_edges:
	ids = line.rstrip().split('.')
	if num_copies > 10:
		id1 = ids[0][-6:]
		id2 = ids[1][:6]
	else:
		id1 = ids[0][-5:]
		id2 = ids[1][:5]

	if id1 in inputs_dict and id2 in outputs_dict:
		arcs_ij.append((inputs_dict[id1],outputs_dict[id2]))

	if id1 in inputs_dict and id2 in pools_dict:
		arcs_il.append((inputs_dict[id1],pools_dict[id2]))

	if id1 in pools_dict and id2 in outputs_dict:
		arcs_lj.append((pools_dict[id1],outputs_dict[id2]))
		

# next, extract the nominal input qualities
# ignore limits on output qualities
specs_in_nominal = [[0.0 for k in range(num_specs)] for i in range(num_inputs)]

for line in lines_specs:
	if num_specs > 1:
		attr = int(line.rstrip().split('.')[0][-1])
	else:
		attr = 0

	if (num_specs == 1) or (attr > 0):
		if num_copies > 10:
			id = line.rstrip().split('.')[1][:6]
			val = float(line[14:])
		else:
			id = line.rstrip().split('.')[1][:5]
			val = float(line[13:])
	else:
		if num_copies > 10:
			id = line.rstrip().split('.')[1][:6]
			val = float(line[16:])
		else:
			id = line.rstrip().split('.')[1][:5]
			val = float(line[15:])

	specs_in_nominal[inputs_dict[id]][attr] = val


# define aliases for problem dimensions
num_i = num_inputs
num_l = num_pools
num_j = num_outputs
num_k = num_specs


# generate auxiliary sets for the model
inputs_to_pool = [[t[0] for t in arcs_il if t[1] == l] for l in range(num_l)]  # which inputs are connected to each pool?
inputs_to_output = [[t[0] for t in arcs_ij if t[1] == j] for j in range(num_j)]  # which inputs are directly connected to each output?
outputs_from_input = [[t[1] for t in arcs_ij if t[0] == i] for i in range(num_i)]  # which outputs are directly connected to each input?
outputs_from_pool = [[t[1] for t in arcs_lj if t[0] == l] for l in range(num_l)]  # which outputs are connected to each pool?
pools_from_input = [[t[1] for t in arcs_il if t[0] == i] for i in range(num_i)] # which pools are connected to each input?
pools_to_output = [[t[0] for t in arcs_lj if t[1] == j] for j in range(num_j)]  # which pools are directly connected to each output?
inputs_pools_to_output = [[(i,t[0]) for t in arcs_lj if t[1] == j for i in inputs_to_pool[t[0]]] for j in range(num_j)]  # which input-pool pairs connected to each output?
arcs_ilj = [(i,l,j) for l in range(num_l) for i in inputs_to_pool[l] for j in outputs_from_pool[l]]  # list of arcs from inputs to pools to outputs

outputs_reachable_from_input = [list(set([t[1] for t in arcs_ij if t[0] == i] + flatten([[t[1] for t in arcs_lj if t[0] == l if i in inputs_to_pool[l]] for l in range(num_l)]))) for i in range(num_i)]  # which outputs are reachable starting from each input?
inputs_reachable_from_output = [[i for i in range(num_i) if j in outputs_reachable_from_input[i]] for j in range(num_j)]  # which inputs are reachable starting from each output reversing arc directions?


# randomly generate target qualities at outputs
specs_out_min = [[0.0 for k in range(num_specs)] for j in range(num_outputs)]
specs_out_max = [[0.0 for k in range(num_specs)] for j in range(num_outputs)]

for j in range(num_j):
	for k in range(num_k):
		spec_out_min = np.amin([specs_in_nominal[i][k] for i in inputs_reachable_from_output[j]])
		spec_out_max = np.amax([specs_in_nominal[i][k] for i in inputs_reachable_from_output[j]])

		start_frac = 0.2 + 0.2*random.random()
		stop_frac = 0.6 + 0.2*random.random()

		specs_out_min[j][k] = spec_out_min + start_frac*(spec_out_max - spec_out_min)
		specs_out_max[j][k] = spec_out_min + stop_frac*(spec_out_max - spec_out_min)



for inst in range(numInstances):

	output_file = instances_path + "pooling_c15_e150_q1_" + str(inst+1) + ".jl"
	json_file = instances_path + "pooling_c15_e150_q1_" + str(inst+1) + ".json"

	theta_specs = np.random.rand(num_inputs, num_specs)
	specs_in = [[0.0 for k in range(num_specs)] for i in range(num_inputs)]
	for i in range(num_inputs):
		for k in range(num_specs):
			specs_in[i][k] = (1 - pert_frac + 2*pert_frac*theta_specs[i,k])*specs_in_nominal[i][k]


	# define auxiliary parameters for the quality constraints
	gamma_low = [[[specs_in[i][k]-specs_out_min[j][k] for k in range(num_specs)] for j in range(num_outputs)] for i in range(num_inputs)]
	gamma_up = [[[specs_in[i][k]-specs_out_max[j][k] for k in range(num_specs)] for j in range(num_outputs)] for i in range(num_inputs)]


	# Write to file. Convert to one indexing before creating Julia input
	f = open(output_file, "w")
	f.write("using SparseArrays \n\n")

	f.write("num_i = " + str(num_i) + "\n")
	f.write("num_l = " + str(num_l) + "\n")
	f.write("num_j = " + str(num_j) + "\n")
	f.write("num_k = " + str(num_k) + "\n\n")

	f.write("arcs_ij = " + str([(i+1,j+1) for (i,j) in arcs_ij]) + "\n")
	f.write("arcs_il = " + str([(i+1,l+1) for (i,l) in arcs_il]) + "\n")
	f.write("arcs_lj = " + str([(l+1,j+1) for (l,j) in arcs_lj]) + "\n")
	f.write("arcs_ilj = " + str([(i+1,l+1,j+1) for (i,l,j) in arcs_ilj]) + "\n\n")

	f.write("inputs_to_pool = [Int64[] for l=1:num_l] \n")
	for l in range(num_l):
		f.write("inputs_to_pool[" + str(l+1) + "] = " + str([i+1 for i in inputs_to_pool[l]]) + "\n")
	f.write("\n")

	f.write("inputs_to_output = [Int64[] for j=1:num_j] \n")
	for j in range(num_j):
		f.write("inputs_to_output[" + str(j+1) + "] = " + str([i+1 for i in inputs_to_output[j]]) + "\n")
	f.write("\n")

	f.write("pools_to_output = [Int64[] for j=1:num_j] \n")
	for j in range(num_j):
		f.write("pools_to_output[" + str(j+1) + "] = " + str([l+1 for l in pools_to_output[j]]) + "\n")
	f.write("\n")

	f.write("pools_from_input = [Int64[] for i=1:num_i] \n")
	for i in range(num_i):
		f.write("pools_from_input[" + str(i+1) + "] = " + str([l+1 for l in pools_from_input[i]]) + "\n")
	f.write("\n")

	f.write("outputs_from_pool = [Int64[] for l=1:num_l] \n")
	for l in range(num_l):
		f.write("outputs_from_pool[" + str(l+1) + "] = " + str([j+1 for j in outputs_from_pool[l]]) + "\n")
	f.write("\n")

	f.write("outputs_from_input = [Int64[] for i=1:num_i] \n")
	for i in range(num_i):
		f.write("outputs_from_input[" + str(i+1) + "] = " + str([j+1 for j in outputs_from_input[i]]) + "\n")
	f.write("\n")

	f.write("inputs_pools_to_output = [Tuple{Int64, Int64}[] for j=1:num_j] \n")
	for j in range(num_j):
		f.write("inputs_pools_to_output[" + str(j+1) + "] = " + str([(i+1,l+1) for (i,l) in inputs_pools_to_output[j]]) + "\n")
	f.write("\n")


	f.write("cap_i = " + str(cap_i) + "\n")
	f.write("cap_l = " + str(cap_l) + "\n")
	f.write("cap_j = " + str(cap_j) + "\n\n")

	f.write("cost_ij_rows = " + str([i+1 for i in cost_ij_rows]) + "\n")
	f.write("cost_ij_cols = " + str([j+1 for j in cost_ij_cols]) + "\n")
	f.write("cost_ij_vals = " + str(cost_ij_vals) + "\n")
	f.write("cost_ij = sparse(cost_ij_rows, cost_ij_cols, cost_ij_vals) \n\n")

	f.write("cost_il_rows = " + str([i+1 for i in cost_il_rows]) + "\n")
	f.write("cost_il_cols = " + str([l+1 for l in cost_il_cols]) + "\n")
	f.write("cost_il_vals = " + str(cost_il_vals) + "\n")
	f.write("cost_il = sparse(cost_il_rows, cost_il_cols, cost_il_vals) \n\n")

	f.write("cost_lj_rows = " + str([l+1 for l in cost_lj_rows]) + "\n")
	f.write("cost_lj_cols = " + str([j+1 for j in cost_lj_cols]) + "\n")
	f.write("cost_lj_vals = " + str(cost_lj_vals) + "\n")
	f.write("cost_lj = sparse(cost_lj_rows, cost_lj_cols, cost_lj_vals) \n\n")

	f.write("gamma_low = " + str(gamma_low) + "\n\n")
	f.write("gamma_up = " + str(gamma_up) + "\n\n\n\n")


	f.write("# ----- Variables ----- #\n")
	f.write("@variables(m, begin \n")
	f.write("	0 <= w[(i,l,j) in arcs_ilj] <= min(cap_i[i],cap_l[l],cap_j[j]) \n")
	f.write("	0 <= x_il[(i,l) in arcs_il] <= min(cap_i[i],cap_l[l]) \n")
	f.write("	0 <= x_lj[(l,j) in arcs_lj] <= min(cap_l[l],cap_j[j]) \n")
	f.write("	0 <= x_ij[(i,j) in arcs_ij] <= min(cap_i[i],cap_j[j]) \n")
	f.write("	0 <= q[arcs_il] <= 1 \n")
	f.write("end) \n\n\n")


	f.write("# ----- Constraints ----- #\n")
	f.write("@constraints(m, begin \n")
	f.write("	capi[i=1:num_i], sum(x_ij[(i,j)] for j in outputs_from_input[i]) + sum(x_il[(i,l)] for l in pools_from_input[i]) <= cap_i[i] \n\n")
	f.write("	capl[l=1:num_l], sum(x_lj[(l,j)] for j in outputs_from_pool[l]) <= cap_l[l] \n\n")
	f.write("	capj[j=1:num_j], sum(x_ij[(i,j)] for i in inputs_to_output[j]) + sum(x_lj[(l,j)] for l in pools_to_output[j]) <= cap_j[j] \n\n")
	f.write("	sumfrac[l=1:num_l], sum(q[(i,l)] for i in inputs_to_pool[l]) == 1 \n\n")
	f.write("	eq1[(i,l) in arcs_il], x_il[(i,l)] == sum(w[(i,l,j)] for j in outputs_from_pool[l]) \n\n")
	f.write("    specdown[k=1:num_k, j=1:num_j], sum(gamma_low[i][j][k]*w[(i,l,j)] for (i,l) in inputs_pools_to_output[j]) + sum(gamma_low[i][j][k]*x_ij[(i,j)] for i in inputs_to_output[j]) >= 0 \n\n")
	f.write("	specup[k=1:num_k, j=1:num_j], sum(gamma_up[i][j][k]*w[(i,l,j)] for (i,l) in inputs_pools_to_output[j]) + sum(gamma_up[i][j][k]*x_ij[(i,j)] for i in inputs_to_output[j]) <= 0 \n\n")
	f.write("	rlt1[(l,j) in arcs_lj], sum(w[(i,l,j)] for i in inputs_to_pool[l]) == x_lj[(l,j)] \n\n")
	f.write("	rlt2[(i,l) in arcs_il], sum(w[(i,l,j)] for j in outputs_from_pool[l]) <= cap_l[l]*q[(i,l)] \n\n")
	f.write("end) \n\n\n")

	f.write("@NLconstraint(m, blin[(i,l,j) in arcs_ilj], w[(i,l,j)] == q[(i,l)]*x_lj[(l,j)]) \n\n\n")


	f.write("# ----- Objective ----- #\n")
	f.write("@objective(m, Min, sum(cost_ij[i,j]*x_ij[(i,j)] for (i,j) in arcs_ij) + sum(cost_il[i,l]*x_il[(i,l)] for (i,l) in arcs_il) + sum(cost_lj[l,j]*x_lj[(l,j)] for (l,j) in arcs_lj) ) \n")
	
		
	f.close()


	json_dict = {}
	json_dict["num_inputs"] = num_i
	json_dict["num_pools"] = num_l
	json_dict["num_outputs"] = num_j
	json_dict["num_qualities"] = num_k

	json_dict["num_copies"] = num_copies
	json_dict["num_added_edges"] = num_added_edges
	json_dict["pert_frac"] = pert_frac
	json_dict["capacities_scaling"] = capacities_scaling
	json_dict["costs_scaling"] = costs_scaling

	json_dict["arcs_inputs_outputs"] = [(i+1,j+1) for (i,j) in arcs_ij]
	json_dict["arcs_inputs_pools"] = [(i+1,l+1) for (i,l) in arcs_il]
	json_dict["arcs_pools_outputs"] = [(l+1,j+1) for (l,j) in arcs_lj]
	json_dict["arcs_inputs_pools_outputs"] = [(i+1,l+1,j+1) for (i,l,j) in arcs_ilj]

	json_dict["capacities_inputs"] = cap_i
	json_dict["capacities_pools"] = cap_l
	json_dict["capacities_outputs"] = cap_j

	json_dict["cost_ij_rows"] = [i+1 for i in cost_ij_rows]
	json_dict["cost_ij_cols"] = [i+1 for i in cost_ij_rows]
	json_dict["cost_ij_vals"] = cost_ij_vals
	json_dict["cost_il_rows"] = [i+1 for i in cost_il_rows]
	json_dict["cost_il_cols"] = [i+1 for i in cost_il_rows]
	json_dict["cost_il_vals"] = cost_il_vals
	json_dict["cost_lj_rows"] = [i+1 for i in cost_lj_rows]
	json_dict["cost_lj_cols"] = [i+1 for i in cost_lj_rows]
	json_dict["cost_lj_vals"] = cost_lj_vals

	json_dict["specs_inputs_nominal"] = specs_in_nominal
	json_dict["theta_specs_inputs"] = theta_specs.tolist()
	json_dict["specs_inputs"] = specs_in
	json_dict["specs_outputs_min"] = spec_out_min
	json_dict["specs_outputs_max"] = spec_out_max
	json_dict["gamma_low"] = gamma_low
	json_dict["gamma_up"] = gamma_up


	with open(json_file, 'w', encoding='utf-8') as outF:
		json.dump(json_dict, outF, ensure_ascii=False, indent=4)
