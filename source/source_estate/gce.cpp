#include "gce.h"

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
    workfunction_ = mixed_workfunction;
    has_workfunction_ = true;
}

} // namespace elecstate
