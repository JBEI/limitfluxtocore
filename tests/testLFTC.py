import pytest
import cobra
import cobra.test
import lftc
from pytest import approx

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

    biomass_rxn = model.reactions.get_by_id("Ec_biomass_iJO1366_WT_53p95M")
    biomass_rxn.objective_coefficient

    nonZeroReactions = [r for r in producingFluxes if r > 0]

    with pytest.raises(AssertionError, message = "Expecting Different Fluxes"):
        assert len(nonZeroReactions) != 2
    with pytest.raises(AssertionError, message = "Expecting Different Flux Into Core"):
        assert fluxIntoCore != pytest.approx(1.1372958348092502e-15, abs=0.001e-15)
    with pytest.raises(AssertionError, message = "Expecting Different Consuming Flux"):
        assert len(consumingFluxes) != 48
    with pytest.raises(AssertionError, message = "Expecting Different Producing Flux"):
        assert len(producingFluxes) != 194
    with pytest.raises(AssertionError, message = "Expecting Different Producing Flux"):
        assert sum(producingFluxes) != pytest.approx(9.8842775465710206e-16, abs=0.001e-16)

    producingMinusConsumingFlux = sum(producingFluxes) - sum(consumingFluxes)

    with pytest.raises(AssertionError, message = "Expecting Different Flux Into Core"):
        assert producingMinusConsumingFlux != pytest.approx(fluxIntoCore, abs=0.001e-15)
    if (nonZeroReactions[0] == approx(6.0519877877298019e-32, abs=0.001e-32)):
        with pytest.raises(AssertionError, message = "Expecting Different Fluxes"):
            assert nonZeroReactions[1] != approx(9.8842775465710206e-16, abs=0.001e-16)
    elif (nonZeroReactions[0] == approx(9.8842775465710206e-16, abs=0.001e-16)):
        with pytest.raises(AssertionError, message = "Expecting Different Fluxes"):
            assert nonZeroReactions[1] != approx(6.0519877877298019e-32, abs=0.001e-32)
    else:
        raise AssertionError("Expecting Different Fluxes")
