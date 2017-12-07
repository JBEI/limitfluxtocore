"""
These are the user level functions and objects for lftc

By Tyler W. H. Backman
"""

import pandas as pd
import numpy as np
import cobra
import re
from .anneal import Annealer
from .exploreModel import \
    findSubsetConnectedToFeed, \
    findProducingReactions, \
    findConsumingReversibleReactions, \
    findConsumedMetabolites, \
    findProducingReactionsOutsideCore, \
    findReversibleConsumingReactionsOutsideCore

# define a default set of currency metabolites
currencyMetabolites = set(['h2o_c', 'pi_c', 'co2_c', 'atp_c', 'coa_c', 'imp_c',
                           'gmp_c', 'adp_c', 'amp_c', 'h_c', 'h_e', 'nadph_c', 
                           'q6h2_c', 'q6_c', 'nadp_c', 'nadh_c', 'nad_c', 
                           'nh4_c', 'ppi_c', 'gtp_c', 'gdp_c', 'o2_c', 
                           'q8h2_c', 'q8_c', 'fad_c', 'fadh2_c'])
currencyMetabolites.update({re.sub('_c', '_m', r) for r in currencyMetabolites})

def limitFluxToCore(
    coreReactionNames, 
    model, 
    currencyMetabolites=currencyMetabolites, 
    copyModel=True
    ):
    """Main limit flux to core algorithm.

    Minimizes the sum of fluxes into core metabolism via linear optimization
    subject to the following options:

    Args:
        coreReactionNames (set): The set of reaction names of type str from 
            model included in core.
        model (cobra.core.model.Model): A COBRApy genome scale model.
        currencyMetabolites (set): Optional, a set of metabolites to exclude 
            when identifying reactions which feed carbon into the core. If 
            excluded, the lftc.currencyMetabolites default set is used.
        copyModel (bool): Should the model be copied first to avoid 
            modification? Otherwise, it's objective function will be altered.

    Returns:
        (tuple): tuple containing:
            arg1 (numpy.float64): The sum of fluxes into core metabolism.
            arg2 (pandas.core.series.Series): The (positive or zero) upper 
                flux bound of all fluxes that produce metabolites in the core.
            arg3 (pandas.core.series.Series): The (negative or zero) lower 
                flux bound of all reversible reaction fluxes that can produce 
                metabolites in the core in reverse direction.
    """

    # sanity check inputs
    assert type(model) is cobra.core.model.Model, \
        'model is not of type cobra.core.model.Model'
    assert type(coreReactionNames) is set, 'coreReactionNames is not a set'
    assert 0 < len(coreReactionNames) <= len(model.reactions), \
        'invalid size for coreReactionNames'
    assert len(coreReactionNames.intersection([r.id for r in model.reactions])) \
        == len(coreReactionNames), 'some core reaction names missing from model'
    assert type(copyModel) is bool, 'copyModel not type bool'
    assert type(currencyMetabolites) is set, 'currencyMetabolites is not a set'

    if copyModel:
        model = model.copy()
    
    # find reactions to minimize
    allMetaboliteNames = findConsumedMetabolites(
        coreReactionNames, 
        model, 
        currencyMetabolites,
        )
    allProducingReactions = findProducingReactionsOutsideCore(
        allMetaboliteNames, 
        model, 
        coreReactionNames,
        )
    allConsumingReactions = findReversibleConsumingReactionsOutsideCore(
        allMetaboliteNames, 
        model, 
        coreReactionNames,
        )
    
    # set the objective function to simultaneously minimize reactions 
    # which produce metabolites in core, and maximize
    # reversible reactions which consume metabolites from core
    objective = {}
    for reactionName in allProducingReactions:
        reaction = model.reactions.get_by_id(reactionName)
        objective[reaction.forward_variable] = 1
    for reactionName in allConsumingReactions:
        reaction = model.reactions.get_by_id(reactionName)
        objective[reaction.reverse_variable] = 1
        
    # set the objective for the model, and run FBA
    model.objective = {}
    model.objective.set_linear_coefficients(objective)
    model.objective.direction = 'min'
    fba_solution = model.optimize()
    
    # get fluxes from solution
    fluxIntoCore = 0
    producingFluxes = fba_solution.fluxes[allProducingReactions]
    for producingFlux in producingFluxes.index:
        if producingFluxes[producingFlux] < 0:
            producingFluxes[producingFlux] = 0
        fluxIntoCore += producingFluxes[producingFlux]
    consumingFluxes = fba_solution.fluxes[allConsumingReactions]
    for consumingFlux in consumingFluxes.index:
        if consumingFluxes[consumingFlux] > 0:
            consumingFluxes[consumingFlux] = 0
        fluxIntoCore -= consumingFluxes[consumingFlux]
    
    return fluxIntoCore, producingFluxes, consumingFluxes

def setModelFluxes(model, producingFluxes, consumingFluxes):
    """Apply fluxes to genome scale model.

    Applies flux limits returned by limitFluxToCore() to a genome scale model:

    Args:
        model (cobra.core.model.Model): A COBRApy genome scale model.
        producingFluxes (pandas.core.series.Series): The (positive or zero) 
            upper flux bound of all fluxes that produce metabolites in the core.
        consumingFluxes (pandas.core.series.Series): The (negative or zero) 
            lower flux bound of all reversible reaction 
            fluxes that can produce metabolites in the core in reverse direction.

    Returns:
        model (cobra.core.model.Model): A COBRApy genome scale model with 
            new flux bounds.
    """

    # sanity check inputs
    assert type(model) is cobra.core.model.Model, \
        'model is not of type cobra.core.model.Model'
    assert type(producingFluxes) is pd.core.series.Series
    assert type(consumingFluxes) is pd.core.series.Series
    reactionNames = [r.id for r in model.reactions]
    assert len(set(producingFluxes.index).intersection(reactionNames)) \
        == len(producingFluxes)
    assert len(set(consumingFluxes.index).intersection(reactionNames)) \
        == len(consumingFluxes)

    newModel = model.copy()

    for reactionName, fluxLimit in producingFluxes.iteritems():
        reaction = newModel.reactions.get_by_id(reactionName)
        reaction.upper_bound = fluxLimit

    for reactionName, fluxLimit in consumingFluxes.iteritems():
        reaction = newModel.reactions.get_by_id(reactionName)
        reaction.lower_bound = fluxLimit

    return newModel

# specify the optimization problem
class OptimalCoreProblem(Annealer):

    def __init__(
        self, 
        state, 
        model,
        feed, 
        currencyMetabolites=currencyMetabolites, 
        minOverlapWithStart=0.8, 
        maxOverlapWithModel=1.0,
        excludeReactions=set(),
        logFile=None,
        ):
        """Simulated Annealing Core Optimizer.

        Applies simulated annealing to iteratively add or remove reactions
        from core, to minimize the required sum of fluxes into core metabolism.
        See docs at https://github.com/perrygeo/simanneal 
        for use of this object once instantiated.

        Args:
            state (set): The set of reaction names of type str to set the 
                initial starting core state.
            model (cobra.core.model.Model): A COBRApy genome scale model.
            feed (str): The carbon uptake feed.
            currencyMetabolites (set): Optional, a set of metabolites to 
                exclude when identifying reactions which feed carbon into the
                core. If excluded, the lftc.currencyMetabolites default
                set is used.
            minOverlapWithStart (float): A float between 0 and 1 indicating 
                the minimum fraction of reactions from the initial state that
                we should require stay in the final solution. This sets
                a lower bound on core size.
            maxOverlapWithModel (float): A float between 0 and 1 indicating 
                the maximum fraction of reactions from the genome scale model
                to include in the core. This sets an upper bound on core size.
            excludeReactions (set): Optional, a set of reaction names of type
                str to exclude from any possible core solutions. Sometimes it
                is desirable to include exchange fluxes here so they don't get
                added to the core.
        """
        # sanity check inputs
        assert type(model) is cobra.core.model.Model, \
            'model is not of type cobra.core.model.Model'
        reactionNames = [r.id for r in model.reactions]
        assert type(state) is set, 'state is not a set'
        assert 0 < len(state) <= len(model.reactions), \
            'invalid size for state'
        assert len(state.intersection([r.id for r in model.reactions])) \
            == len(state), 'some initial state reaction names missing from model'
        assert type(feed) is str
        assert feed in reactionNames
        assert type(currencyMetabolites) is set
        assert type(minOverlapWithStart) is float
        assert 0 <= minOverlapWithStart <= 1
        assert type(maxOverlapWithModel) is float
        assert 0 <= maxOverlapWithModel <= 1
        assert type(excludeReactions) is set
        assert len(set(excludeReactions).intersection(reactionNames)) \
            == len(excludeReactions)

        
        self.currencyMetabolites = currencyMetabolites
        self.model = model.copy()
        self.maxSize = float(maxOverlapWithModel) * len(self.model.reactions)
        self.feed = feed

        # confirm that initial state is connected
        connected = findSubsetConnectedToFeed(
            state,
            self.feed,
            self.currencyMetabolites,
            self.model,
            )
        if len(connected) < len(state):
            print('warning, initial core not connected, keeping only', 
                len(connected), 'reactions out of', len(state))
            state = connected

        # save initial state of core
        self.startSet = state
        fluxIntoCore, producingFluxes, consumingFluxes = limitFluxToCore(
            self.startSet, 
            self.model, 
            self.currencyMetabolites,
            False,
            ) 
        self.producingFluxes = producingFluxes
        self.consumingFluxes = consumingFluxes

        # initialize logfile
        if logFile:
            self.logFile = open(logFile, 'a')
            self.logFile.write('energy,size\n')
        else:
            self.logFile = None
        
        self.excludeReactions = excludeReactions
        self.minOverlapWithStart = float(minOverlapWithStart)

        assert self.minOverlapWithStart*len(self.startSet) < \
            self.maxSize, 'max size not greater than min size!'
        super(OptimalCoreProblem, self).__init__(state)  # important!

    def move(self):
        # randomly adds or removes a reaction from the core

        # update the logfile with temp from previous move
        if self.logFile:
            self.logFile.write(str(self.energy()) + ',' + \
                str(len(self.state)) + '\n')

        if np.random.randint(2) and (len(self.state) < self.maxSize):
            # half of the time add a reaction
            boundaryReactionNames = set(self.producingFluxes.index)
            boundaryReactionNames.update(self.consumingFluxes.index)
            
            boundaryReactionNames.difference_update(self.excludeReactions)
            newReaction = np.random.choice(tuple(boundaryReactionNames))
            self.state.update([newReaction])
        else:
            # the other half of the time,
            # remove a reaction, and double check that the core remains connected
            stateCopy = self.state.copy()

            # if we hit the minimum overlap with the start, don't remove any more
            # reactions that overlap with start
            currentOverlapLength = len(self.startSet.intersection(self.state))
            overlapWithStart = currentOverlapLength / len(self.startSet)
            if overlapWithStart <= self.minOverlapWithStart:
                stateCopy.difference_update(self.startSet)
                
            # randomize the set of reactions to potentially remove
            stateList = list(stateCopy)
            np.random.shuffle(stateList)
                
            for reactionName in stateList:
                tempState = self.state.copy()
                tempState.remove(reactionName)
                connected = findSubsetConnectedToFeed(
                    tempState,
                    self.feed,
                    self.currencyMetabolites,
                    self.model,
                    )
                if len(connected) == len(tempState):
                    self.state = connected
                    break

    def energy(self):
        # calculate the score

        fluxIntoCore, self.producingFluxes, self.consumingFluxes = \
            limitFluxToCore(
                self.state, 
                self.model, 
                self.currencyMetabolites, 
                False,
                )

        return fluxIntoCore

    def prune(self):
        # Removes any individual newly added reactions which don't
        # reduce total energy into the core. This reverses the addition
        # of unnecessary individual reactions, however it does not
        # combinatorially test removing multiple reactions simultaneously.

        newReactions = self.state.difference(self.startSet)

        for newReaction in newReactions:
            oldState = self.state.copy()

            # try removing a reaction and test if the core is still connected
            tempState = self.state.copy()
            tempState.remove(newReaction)
            connected = findSubsetConnectedToFeed(
                tempState,
                self.feed,
                self.currencyMetabolites,
                self.model,
                )

            # if core is still connected, check the new energy
            if len(connected) == len(tempState):
                currentEnergy = self.energy()
                self.state = connected
                newEnergy = self.energy()

                # if new energy is worse, go back to the old solution
                # otherwise keep it
                if newEnergy > currentEnergy:
                    self.state = oldState
