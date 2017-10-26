"""
These are backend functions for navigating an SBML model using COBRApy

By Tyler W. H. Backman
"""

def findSubsetConnectedToFeed(coreReactionNames, feed, currencyMetabolites, model):
    # start with the feed and make sure core is connected to feed

    newCoreReactionNames = set()
    unexploredReactionNames = {feed} # this starts with just the feed name
    while len(unexploredReactionNames) > 0:
        reactionName = unexploredReactionNames.pop()
        newCoreReactionNames.update([reactionName])
        reaction = model.reactions.get_by_id(reactionName)
        connectedMetabolites = {m.id for m in reaction.reactants}
        connectedMetabolites.update({m.id for m in reaction.products})
        connectedMetabolites = connectedMetabolites.difference(currencyMetabolites)

        for metaboliteName in connectedMetabolites:
            metabolite = model.metabolites.get_by_id(metaboliteName)
            connectedReactionNames = {r.id for r in metabolite.reactions}
            connectedReactionNames = connectedReactionNames.intersection(coreReactionNames)
            newReactionsThisLoop = connectedReactionNames.difference(newCoreReactionNames)
            unexploredReactionNames.update(newReactionsThisLoop)
            newCoreReactionNames.update(newReactionsThisLoop)
    return newCoreReactionNames
        
def findProducingReactions(metabolite):
    # find all reactions which can produce a given metabolite
    
    allReactions = list(metabolite.reactions)
    producingReactions = []
    
    for reaction in allReactions:
        if metabolite in reaction.products:
            producingReactions.append(reaction)
            
    return producingReactions

def findConsumingReversibleReactions(metabolite):
    # find all reactions which can consume a given metabolite and are reversible
    
    allReactions = list(metabolite.reactions)
    allReactions = list(filter(lambda r: r.reversibility, allReactions))
    consumingReactions = []
    
    for reaction in allReactions:
        if metabolite in reaction.reactants:
            consumingReactions.append(reaction)
            
    return consumingReactions

def findConsumedMetabolites(coreReactionNames, model, currencyMetabolites):
    # find all potential sources and exchanges: e.g. metabolites which can be consumed by core reactions

    allMetaboliteNames = set()
    for reactionName in coreReactionNames:
        reaction = model.reactions.get_by_id(reactionName)
        allMetaboliteNames.update({m.id for m in reaction.reactants})
        if reaction.reversibility:
            allMetaboliteNames.update({m.id for m in reaction.products})
        
    # exclude currency reactions that don't transfer any carbon into core
    allMetaboliteNames = allMetaboliteNames.difference(currencyMetabolites)
    return allMetaboliteNames

def findProducingReactionsOutsideCore(allMetaboliteNames, model, coreReactionNames):
    # find all reactions which can produce these metabolites, and are not in the core

    allProducingReactions = set()
    for metaboliteName in list(allMetaboliteNames):
        metabolite = model.metabolites.get_by_id(metaboliteName)
        producingReactions = findProducingReactions(metabolite)
        externalProducingReactions = {r.id for r in producingReactions}.difference(coreReactionNames)
        allProducingReactions.update(externalProducingReactions)
    return allProducingReactions

def findReversibleConsumingReactionsOutsideCore(allMetaboliteNames, model, coreReactionNames):
    # find all reactions which consume metabolites from core, but are reversible, so they can
    # potentially feed carbon into the core

    allConsumingReactions = set()
    for metaboliteName in list(allMetaboliteNames):
        metabolite = model.metabolites.get_by_id(metaboliteName)
        consumingReactions = findConsumingReversibleReactions(metabolite)
        externalConsumingReactions = {r.id for r in consumingReactions}.difference(coreReactionNames)
        allConsumingReactions.update(externalConsumingReactions)
    return allConsumingReactions
