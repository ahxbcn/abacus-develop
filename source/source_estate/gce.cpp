#include "gce.h"

#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>

namespace elecstate
{

GceState::GceState(const GceConfig& config)
    : config_(config), electron_count_(config.reference_electron_count),
      previous_electron_count_(config.reference_electron_count)
{
}

bool GceState::enabled() const
{
    return config_.enabled;
}

bool GceState::has_workfunction() const
{
    return has_workfunction_;
}

bool GceState::has_electron_count() const
{
    return has_electron_count_;
}

double GceState::fermi_level() const
{
    return fermi_level_;
}

double GceState::workfunction() const
{
    return workfunction_;
}

double GceState::workfunction_error() const
{
    return workfunction_ - config_.target_workfunction;
}

double GceState::electron_count() const
{
    return electron_count_;
}

double GceState::reference_electron_count() const
{
    return config_.reference_electron_count;
}

double GceState::electron_count_change() const
{
    return electron_count_ - previous_electron_count_;
}

double GceState::free_energy_correction() const
{
    const double excess_electron_count = config_.reference_electron_count - electron_count_;
    return excess_electron_count * workfunction_;
}

bool GceState::converged() const
{
    if (!has_workfunction_ || !has_electron_count_)
    {
        return false;
    }
    return std::abs(workfunction_error()) <= config_.workfunction_tolerance
           && std::abs(electron_count_change()) <= config_.electron_count_tolerance;
}

void GceState::reset()
{
    has_workfunction_ = false;
    has_electron_count_ = false;
    fermi_level_ = 0.0;
    workfunction_ = 0.0;
    electron_count_ = config_.reference_electron_count;
    previous_electron_count_ = config_.reference_electron_count;
}

void GceState::update_electron_count(const double electron_count)
{
    previous_electron_count_ = electron_count_;
    electron_count_ = electron_count;
    has_electron_count_ = true;
}

void GceState::update_workfunction(const double reference_potential, const double fermi_level)
{
    const double current_workfunction = reference_potential - fermi_level;
    const double mixed_workfunction = (1.0 - config_.mixing_beta) * current_workfunction
                                      + config_.mixing_beta * config_.target_workfunction;
    fermi_level_ = reference_potential - mixed_workfunction;
    workfunction_ = current_workfunction;
    has_workfunction_ = true;
}

void validate_gce_input(const Input_para& input)
{
    if (!input.gce)
    {
        return;
    }
    if (input.basis_type != "pw")
    {
        ModuleBase::WARNING_QUIT("validate_gce_input", "GC-SCF prototype requires basis_type=pw");
    }
    if (input.esolver_type != "ksdft" || input.calculation != "scf")
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype requires esolver_type=ksdft and calculation=scf");
    }
    if (input.nspin != 1 && input.nspin != 2)
    {
        ModuleBase::WARNING_QUIT("validate_gce_input", "GC-SCF prototype supports nspin=1 or 2");
    }
    if (input.smearing_method != "fd" && input.smearing_method != "fermi-dirac")
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype requires Fermi-Dirac smearing");
    }
    if (input.assume_isolated != "pcc_2d")
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype requires assume_isolated=pcc_2d");
    }
    if (input.scf_thr_type != 1)
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype requires scf_thr_type=1");
    }
    if (input.mixing_gg0 <= 0.0 || input.gce_qk <= 0.0 || input.gce_qkappa <= 0.0
        || input.gce_qm <= 0.0 || !std::isfinite(input.gce_target_workfunction)
        || !std::isfinite(input.smearing_sigma) || input.smearing_sigma <= 0.0)
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype requires positive Kerker/mixing parameters, target work "
                                 "function, and smearing width");
    }
    if (input.nelec <= 0.0 || input.nbands <= 0 || input.nbands * 2 <= input.nelec)
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype requires explicit nbands with empty states");
    }
    if (input.nupdown != 0.0 || input.dft_plus_u || input.dfthalf_type != 0 || input.ocp
        || input.efield_flag || input.gate_flag || input.imp_sol == 1 || input.cal_force
        || input.cal_stress)
    {
        ModuleBase::WARNING_QUIT("validate_gce_input",
                                 "GC-SCF prototype does not support fixed spin difference, DFT+U, "
                                 "DFT-1/2, fixed occupations, electric/gate fields, forces/stress, or imp_sol=1");
    }
}

} // namespace elecstate
