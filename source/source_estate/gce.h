#ifndef GCE_H
#define GCE_H

class Charge;
struct Input_para;
class UnitCell;

namespace ModulePW
{
class PW_Basis;
}

namespace elecstate
{

class Potential;

struct GceConfig
{
    bool enabled = false;
    double reference_electron_count = 0.0;
    double target_workfunction = 0.0;
    double mixing_beta = 0.3;
    double workfunction_tolerance = 1.0e-4;
    double electron_count_tolerance = 1.0e-6;
};

/**
 * @brief SCF state for a constant-work-function grand-canonical calculation.
 *
 * Energies and potentials use Ry. The reference electron count is the fixed
 * electron count of the input charge state; the current electron count is
 * allowed to fluctuate around it.
 */
class GceState
{
  public:
    explicit GceState(const GceConfig& config);

    bool enabled() const;
    bool has_workfunction() const;
    bool has_electron_count() const;
    double fermi_level() const;
    double workfunction() const;
    double workfunction_error() const;
    double electron_count() const;
    double reference_electron_count() const;
    double electron_count_change() const;
    double free_energy_correction() const;
    bool converged() const;

    void reset();
    void update_electron_count(double electron_count);
    void update_workfunction(double reference_potential, double fermi_level);

  private:
    GceConfig config_;
    bool has_workfunction_ = false;
    bool has_electron_count_ = false;
    double fermi_level_ = 0.0;
    double workfunction_ = 0.0;
    double electron_count_ = 0.0;
    double previous_electron_count_ = 0.0;
};

/**
 * @brief Return the planar-averaged electrostatic potential at the vacuum
 *        minimum for a slab. All values returned by this function use Ry.
 */
double electrostatic_vacuum_reference(const Charge& charge,
                                      const Potential& potential,
                                      const UnitCell& cell,
                                      int axis);

void validate_gce_input(const Input_para& input);

} // namespace elecstate

#endif
