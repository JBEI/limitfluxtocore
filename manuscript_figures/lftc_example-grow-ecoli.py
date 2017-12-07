#!/usr/bin/python3

# run simulated annealing in parallel to test expansion of the standard E. coli core
# 2017 Tyler W. H. Backman 

import multiprocessing
import lftc
import cobra
import pickle
import re

# Import SBML model
model = cobra.io.read_sbml_model('test_data/EciJR904TKs_sbml3.xml')

# set exchange fluxes to experimental values from Toya2010 wt5h 
model.reactions.EX_glc_e_.lower_bound = -11.7
model.reactions.EX_glc_e_.upper_bound = -11.7
model.reactions.EX_ac_e_.lower_bound = 4.3
model.reactions.EX_ac_e_.upper_bound = 4.3
model.reactions.BiomassEcoli.lower_bound = 0.83
model.reactions.BiomassEcoli.upper_bound = 0.89

# read in E. coli core
with open('test_data/REACTIONSwt5h.txt') as f:
    allLines = f.readlines()
reactionLines = list(filter(lambda l: re.match('^\w', l), allLines))
coreReactionNamesFromFile = set([re.sub('\s.*$', '', l) for l in reactionLines])
coreReactionNamesFromFile = set([re.sub('[\(,\)]', '_', l) for l in coreReactionNamesFromFile])
coreReactionNames = coreReactionNamesFromFile.intersection([r.id for r in model.reactions])

print(str(len(coreReactionNamesFromFile)) + ' core reactions in transitions file')
print(str(len(coreReactionNames)) + ' overlapping with model')
print('missing from model: ', coreReactionNamesFromFile.difference(coreReactionNames))

# run simulated annealing in parallel
def anneal(seed):
    auto_schedule = {'steps': 200000, 'tmax': 50000, 'tmin': 0.01, 'updates': 1000}
    ocp = lftc.OptimalCoreProblem(
        state=coreReactionNames,
        model=model, 
        feed='EX_glc_e_',
        minOverlapWithStart=1.0,
        maxOverlapWithModel=0.13,
	logFile='test_data/logs/' + str(seed) + '.csv',
        excludeReactions=exchanges)
    ocp.set_schedule(auto_schedule)
    coreReactions, score = ocp.anneal(seed=seed)
    ocp.logFile.close()
    return [coreReactions, score]

exchanges = {r.id for r in model.exchanges}
cpuCores = 64
seeds = list(range(1,cpuCores*5,5))
pool = multiprocessing.Pool(processes=cpuCores)
results = pool.map(anneal, seeds)

# print best score and core size
print('best core is', min([r[1] for r in results]))

# save results to file
with open('test_data/results_e_coli.p', 'wb') as f:
    pickle.dump(results, f, pickle.HIGHEST_PROTOCOL)
