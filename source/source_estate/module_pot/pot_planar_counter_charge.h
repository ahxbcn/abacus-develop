#ifndef POT_PLANAR_COUNTER_CHARGE_H
#define POT_PLANAR_COUNTER_CHARGE_H

#include "pot_base.h"

#include <vector>

namespace elecstate
{

/// Gaussian counter-charge sheets on both sides of a slab.
///
/// The total counter charge is chosen to cancel the DFT charge
/// (ions minus electrons). The resulting electrostatic potential is added to
/// the electronic effective potential. The energy correction contains the ion
/// cross term and the counter-charge self energy; the electron cross term is
/// already present in the band energy.
class PotPlanarCounterCharge : public PotBase
{
  public:
    PotPlanarCounterCharge(const ModulePW::PW_Basis* basis,
                           double electron_count,
                           int open_axis,
                           double distance_from_edge,
                           double gaussian_width);

    void cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential) override;
    double get_energy() const override;
    const std::vector<double>& electron_potential() const;
    void set_electron_count(double electron_count);

  private:
    double electron_count_;
    int open_axis_;
    double distance_from_edge_;
    double gaussian_width_;
    std::vector<double> electron_potential_;
    double energy_rydberg_ = 0.0;
};

} // namespace elecstate

#endif
