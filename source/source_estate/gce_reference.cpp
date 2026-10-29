#include "gce.h"

#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/h_hartree_pw.h"
#include "source_estate/module_pot/potential_new.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace elecstate
{

namespace
{

std::vector<double> average_along_axis(const ModulePW::PW_Basis& basis,
                                       const std::vector<double>& values,
                                       const int axis,
                                       const bool absolute)
{
    if (axis < 0 || axis > 2)
    {
        return std::vector<double>();
    }

    int length = basis.nz;
    if (axis == 0)
    {
        length = basis.nx;
    }
    else if (axis == 1)
    {
        length = basis.ny;
    }

    std::vector<double> average(length, 0.0);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        int index = 0;
        if (axis == 0)
        {
            index = ir / (basis.ny * basis.nplane);
        }
        else if (axis == 1)
        {
            const int i = ir / (basis.ny * basis.nplane);
            index = ir / basis.nplane - i * basis.ny;
        }
        else
        {
            index = ir % basis.nplane + basis.startz_current;
        }
        if (index < 0 || index >= length)
        {
            return std::vector<double>();
        }
        const double value = absolute ? std::abs(values[ir]) : values[ir];
        average[index] += value;
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(average.data(), length);
#endif

    const int surface = basis.nxyz / length;
    if (surface <= 0)
    {
        return std::vector<double>();
    }
    for (int i = 0; i < length; ++i)
    {
        average[i] /= static_cast<double>(surface);
    }
    return average;
}

} // namespace

double electrostatic_vacuum_reference(const Charge& charge,
                                      const Potential& potential,
                                      const UnitCell& cell,
                                      const int axis)
{
    const ModulePW::PW_Basis* basis = potential.get_rho_basis();
    if (basis == nullptr || basis->nrxx <= 0 || charge.rho == nullptr || axis < 0 || axis > 2)
    {
        return 0.0;
    }

    std::vector<double> density(basis->nrxx, 0.0);
    for (int is = 0; is < charge.nspin; ++is)
    {
        for (int ir = 0; ir < basis->nrxx; ++ir)
        {
            density[ir] += charge.rho[is][ir];
        }
    }

    const ModuleBase::matrix hartree = H_Hartree_pw::v_hartree(
        cell, const_cast<ModulePW::PW_Basis*>(basis), charge.nspin, charge.rho);
    const double* fixed_potential = potential.get_fixed_v();
    const std::vector<double>* pcc_potential = potential.pcc_electron_potential();
    const std::vector<double>* solvent_potential = potential.find_sol_pot();
    if (fixed_potential == nullptr || pcc_potential == nullptr)
    {
        return 0.0;
    }

    std::vector<double> electrostatic(basis->nrxx, 0.0);
    for (int ir = 0; ir < basis->nrxx; ++ir)
    {
        const double solvent_value = solvent_potential == nullptr ? 0.0 : (*solvent_potential)[ir];
        electrostatic[ir] = fixed_potential[ir] + hartree(0, ir) + (*pcc_potential)[ir] + solvent_value;
    }

    const std::vector<double> density_profile = average_along_axis(*basis, density, axis, true);
    const std::vector<double> potential_profile = average_along_axis(*basis, electrostatic, axis, false);
    if (density_profile.empty() || density_profile.size() != potential_profile.size())
    {
        return 0.0;
    }

    const int length = static_cast<int>(density_profile.size());
    const double windows[7] = {0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1};
    double minimum = 1.0e100;
    int minimum_index = 0;
    for (int i = 0; i < length; ++i)
    {
        double smoothed = 0.0;
        for (int win = 0; win < 7; ++win)
        {
            const int index = (i - 3 + win + length) % length;
            smoothed += windows[win] * density_profile[index];
        }
        if (smoothed < minimum)
        {
            minimum = smoothed;
            minimum_index = i;
        }
    }
    return potential_profile[minimum_index];
}

} // namespace elecstate
