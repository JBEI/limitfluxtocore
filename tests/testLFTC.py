import pytest
import cobra
import cobra.test
from lftc import limitFluxToCore

def test_limitFluxToCore():
    model = cobra.test.create_test_model("ecoli")
    coreModel = cobra.test.create_test_model("textbook")
    coreSet = {r.id for r in coreModel.reactions}
    allSet = {r.id for r in model.reactions}
    exampleCore = allSet.intersection(coreSet)

    exchanges = [r.id for r in model.exchanges]
    model.reactions.EX_glc_e.lower_bound = -11.7
    model.reactions.EX_glc_e.upper_bound = 11.7

    fluxIntoCore, producingFluxes, consumingFluxes = lftc.limitFluxToCore(exampleCore, model)

    nonZeroReactions = [r for r in producingFluxes if r > 0]
    assert(len(nonZeroReactions)) == 2
    with pytest.raises(AssertionError):
        len(nonZeroReactions) != 2
    assert nonZeroReactions[0] == pytest.approx(6.0519877877298019e-32, abs=0.001e-32), 'Different Reaction Fluxes'
    with pytest.raises(AssertionError):
        nonZeroReactions[0]!= pytest.approx(6.0519877877298019e-32, abs=0.001e-32)
    assert nonZeroReactions[1] == pytest.approx(9.8842775465710206e-16, abs=0.001e-16),'Different Reaction Fluxes'
    with pytest.raises(AssertionError):
        nonZeroReactions[1]!= pytest.approx(9.8842775465710206e-16, abs=0.001e-16)
    assert fluxIntoCore == pytest.approx(1.1372958348092502e-15, abs=0.001e-15), 'Different fluxIntoCore'
    with pytest.raises(AssertionError):
        fluxIntoCore != pytest.approx(1.1372958348092502e-15, abs=0.001e-15)
    assert len(consumingFluxes) == 48, 'Different consumingFlux'
    with pytest.raises(AssertionError):
        len(consumingFluxes) != 48
    assert len(producingFLuxes) == 194, 'Different producingFlux'
    with pytest.raises(AssertionError):
        len(producingFluxes) != 194
    assert sum(producingFluxes) == pytest.approx(9.8842775465710206e-16, abs=0.001e-16), 'Different producingFlux'
    with pytest.raises(AssertionError):
        sum(producingFluxes) != pytest.approx(9.8842775465710206e-16, abs=0.001e-16)
    producingMinusConsumingFlux = sum(producingFluxes) - sum(consumingFluxes)
    assert producingMinusConsumingFlux == pytest.approx(fluxIntoCore, abs=0.001e-15), 'Incorrect FluxIntoCore'
    with pytest.raises(AssertionError):
        producingMinusConsumingFlux != pytest.approx(fluxIntoCore, abs=0.001e-15)
