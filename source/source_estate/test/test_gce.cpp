#include "gtest/gtest.h"

#include "source_estate/gce.h"
#include "source_estate/occupy.h"

#include <cmath>
#include <vector>

TEST(GceStateTest, WorkfunctionMixingAndElectronBookkeeping)
{
    elecstate::GceConfig config;
    config.enabled = true;
    config.reference_electron_count = 10.0;
    config.target_workfunction = 4.0;
    config.mixing_beta = 0.5;
    config.workfunction_tolerance = 1.0e-4;
    config.electron_count_tolerance = 1.0e-6;

    elecstate::GceState state(config);
    EXPECT_TRUE(state.enabled());
    EXPECT_FALSE(state.has_workfunction());
    EXPECT_FALSE(state.converged());

    state.update_electron_count(9.8);
    state.update_workfunction(-2.0, -7.0);
    EXPECT_TRUE(state.has_workfunction());
    EXPECT_NEAR(state.workfunction(), 4.5, 1.0e-14);
    EXPECT_NEAR(state.fermi_level(), -6.5, 1.0e-14);
    EXPECT_NEAR(state.electron_count_change(), -0.2, 1.0e-14);
    EXPECT_NEAR(state.free_energy_correction(), 0.9, 1.0e-14);

    state.update_electron_count(9.8);
    state.update_workfunction(-2.2, -6.7);
    EXPECT_NEAR(state.workfunction(), 4.25, 1.0e-14);
    EXPECT_NEAR(state.fermi_level(), -6.45, 1.0e-14);
    EXPECT_FALSE(state.converged());
}

TEST(GceStateTest, ResetKeepsImmutableConfiguration)
{
    elecstate::GceConfig config;
    config.enabled = true;
    config.reference_electron_count = 3.0;
    config.target_workfunction = 2.0;
    elecstate::GceState state(config);

    state.update_electron_count(2.5);
    state.update_workfunction(1.0, -2.0);
    state.reset();
    EXPECT_FALSE(state.has_workfunction());
    EXPECT_DOUBLE_EQ(state.electron_count(), 3.0);
    EXPECT_DOUBLE_EQ(state.workfunction(), 0.0);
}

TEST(GceOccupationTest, FixedChemicalPotentialFermiDiracWeights)
{
    ModuleBase::matrix ekb(1, 2);
    ekb(0, 0) = -1.0;
    ekb(0, 1) = 0.0;
    ModuleBase::matrix wg(1, 2);
    const std::vector<double> wk = {2.0};
    const std::vector<int> isk = {0};

    double nelec = 0.0;
    double demet = 0.0;
    Occupy::fweights(1, wk, 2, -0.5, 0.1, -99, ekb, wg, -1, isk, nelec, demet);

    EXPECT_NEAR(wg(0, 0), 1.9866142981514303, 1.0e-14);
    EXPECT_NEAR(wg(0, 1), 0.0133857018485697, 1.0e-14);
    EXPECT_NEAR(nelec, 2.0, 1.0e-14);
    EXPECT_TRUE(std::isfinite(demet));
}
