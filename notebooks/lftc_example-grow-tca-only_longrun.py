#!/usr/bin/python3

# run simulated annealing in parallel to test expansion of TCA cycle

import multiprocessing
import lftc
import cobra
import pickle

# Import SBML model
model = cobra.io.read_sbml_model('test_data/EciJR904TKs_sbml3.xml')

# set exchange fluxes to experimental values from Toya2010 wt5h 
model.reactions.EX_glc_e_.lower_bound = -11.7
model.reactions.EX_glc_e_.upper_bound = -11.7
model.reactions.EX_ac_e_.lower_bound = 4.3
model.reactions.EX_ac_e_.upper_bound = 4.3
model.reactions.BiomassEcoli.lower_bound = 0.83
model.reactions.BiomassEcoli.upper_bound = 0.89

# set initial core to just include glycolysis and TCA cycle
coreReactionNames = {
    'GLCpts',
    'PGI',
    'PFK',
    'FBA',
    'GAPD',
    'PGK',
    'PGM',
    'ENO',
    'PYK',
    'DHAPT',
    'PDH',
    'CS',
    'ACONT',
    'ICL',
    'ICDHyr',
    'AKGDH',
    'SUCOAS',
    'SUCD1i',
    'FUM',
    'MDH',
}

# run simulated annealing in parallel
def anneal(seed):
    auto_schedule = {'steps': 200000, 'tmax': 50000, 'tmin': 0.01, 'updates': 1000}
    ocp = lftc.OptimalCoreProblem(
        state=coreReactionNames,
        model=model, 
        feed='EX_glc_e_',
        minOverlapWithStart=1.0,
        maxOverlapWithModel=0.13,
	logFile='logs/' + str(seed) + '.csv',
        excludeReactions=exchanges)
    ocp.set_schedule(auto_schedule)
    coreReactions, score = ocp.anneal(seed=seed)
    ocp.logFile.close()
    return [coreReactions, score]

exchanges = [r.id for r in model.exchanges]
cpuCores = 64
seeds = list(range(1,cpuCores*5,5))
pool = multiprocessing.Pool(processes=cpuCores)
results = pool.map(anneal, seeds)

# save results to file
with open('results_tca_only_longrun.p', 'wb') as f:
    pickle.dump(results, f, pickle.HIGHEST_PROTOCOL)
