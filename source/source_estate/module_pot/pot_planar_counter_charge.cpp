#include "pot_planar_counter_charge.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_cell/cell_geometry.h"
#include "source_cell/cell_tools.h"
#include "source_estate/module_charge/charge.h"

#include <cmath>
#include <complex>
#include <vector>

namespace elecstate
{

PotPlanarCounterCharge::PotPlanarCounterCharge(const ModulePW::PW_Basis* basis,
                                               const double electron_count,
                                               const int open_axis,
                                               const double distance_from_edge,
                                               const double gaussian_width)
    : electron_count_(electron_count), open_axis_(open_axis),
      distance_from_edge_(distance_from_edge), gaussian_width_(gaussian_width)
{
    this->rho_basis_ = basis;
    this->dynamic_mode = true;
}

void PotPlanarCounterCharge::set_electron_count(const double electron_count)
{
    electron_count_ = electron_count;
}

double PotPlanarCounterCharge::get_energy() const
{
    return energy_rydberg_;
}

const std::vector<double>& PotPlanarCounterCharge::electron_potential() const
{
    return electron_potential_;
}

void PotPlanarCounterCharge::cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential)
{
    ModuleBase::timer::start("PotPlanarCounterCharge", "cal_v_eff");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell->atoms, cell->ntype, cell->lat0);
    const double tolerance = 1.0e-6;
    unitcell::SlabCell slab;
    const bool perpendicular = unitcell::make_slab_cell(cell->latvec, cell->lat0, open_axis_, tolerance, slab);
    if (!perpendicular)
    {
        ModuleBase::WARNING_QUIT("PotPlanarCounterCharge",
                                 "Planar counter charge requires a perpendicular open lattice vector");
    }
    double ionic_charge = 0.0;
    ion_charges_.clear();
    ion_gradients_.clear();
    for (const unitcell::AtomData& atom : atoms)
    {
        ionic_charge += atom.valence_charge;
        ion_charges_.push_back(atom.valence_charge);
        ion_gradients_.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    }
    const double counter_charge = electron_count_ - ionic_charge;
    if (std::abs(counter_charge) < 1.0e-14)
    {
        electron_potential_.assign(basis.nrxx, 0.0);
        energy_rydberg_ = 0.0;
        ModuleBase::timer::end("PotPlanarCounterCharge", "cal_v_eff");
        return;
    }

    const double width = gaussian_width_;
    const double half_length = 0.5 * slab.length;
    if (width <= 0.0 || distance_from_edge_ <= 0.0 || distance_from_edge_ >= half_length)
    {
        ModuleBase::WARNING_QUIT("PotPlanarCounterCharge",
                                 "Gaussian planar counter charge requires a positive width and a plane inside the cell");
    }
    const double center_offset = half_length - distance_from_edge_;
    const double normalization = counter_charge / (std::sqrt(ModuleBase::TWO_PI) * width * slab.area);
    std::vector<ModuleBase::Vector3<double>> grid_positions;
    ModulePW::grid_positions(basis, cell->latvec, cell->lat0, grid_positions);
    std::vector<double> counter_density(basis.nrxx, 0.0);
    double integrated_charge = 0.0;
    const double dv = cell->omega / basis.nxyz;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double coordinate = grid_positions[ir] * slab.normal - slab.origin;
        const double upper = coordinate - center_offset;
        const double lower = coordinate + center_offset;
        const double gaussian = 0.5 * (std::exp(-upper * upper / (2.0 * width * width))
                                       + std::exp(-lower * lower / (2.0 * width * width)));
        counter_density[ir] = normalization * gaussian;
        integrated_charge += counter_density[ir] * dv;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(integrated_charge);
#endif
    if (std::abs(integrated_charge) > 1.0e-14)
    {
        const double scale = counter_charge / integrated_charge;
        for (double& value : counter_density)
        {
            value *= scale;
        }
    }

    std::vector<std::complex<double>> density_g(basis.nmaxgr, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> potential_g(basis.npw, std::complex<double>(0.0, 0.0));
    basis.real2recip(counter_density.data(), density_g.data());
    double self_energy = 0.0;
    const double e2_four_pi = ModuleBase::e2 * ModuleBase::FOUR_PI;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (ig == basis.ig_gge0)
        {
            continue;
        }
        const double gg = basis.gg[ig];
        const double factor = e2_four_pi / (cell->tpiba2 * gg);
        potential_g[ig] = density_g[ig] * factor;
        self_energy += 0.5 * cell->omega * std::norm(density_g[ig]) * factor;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(self_energy);
#endif

    double ion_cross_energy = 0.0;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (ig == basis.ig_gge0)
        {
            continue;
        }
        for (std::size_t iatom = 0; iatom < atoms.size(); ++iatom)
        {
            const unitcell::AtomData& atom = atoms[iatom];
            const double phase = basis.gcar[ig] * atom.position;
            const std::complex<double> phase_factor(std::cos(phase), std::sin(phase));
            const std::complex<double> value = potential_g[ig] * phase_factor;
            ion_cross_energy += atom.valence_charge * value.real();
            ion_gradients_[iatom] -= basis.gcar[ig] * value.imag();
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(ion_cross_energy);
    if (!ion_gradients_.empty())
    {
        std::vector<double> gradient_buffer(ion_gradients_.size() * 3, 0.0);
        for (std::size_t iatom = 0; iatom < ion_gradients_.size(); ++iatom)
        {
            gradient_buffer[3 * iatom] = ion_gradients_[iatom].x;
            gradient_buffer[3 * iatom + 1] = ion_gradients_[iatom].y;
            gradient_buffer[3 * iatom + 2] = ion_gradients_[iatom].z;
        }
        Parallel_Reduce::reduce_pool(gradient_buffer.data(), static_cast<int>(gradient_buffer.size()));
        for (std::size_t iatom = 0; iatom < ion_gradients_.size(); ++iatom)
        {
            ion_gradients_[iatom] = ModuleBase::Vector3<double>(
                gradient_buffer[3 * iatom], gradient_buffer[3 * iatom + 1], gradient_buffer[3 * iatom + 2]);
        }
    }
#endif

    electron_potential_.assign(basis.nrxx, 0.0);
    basis.recip2real(potential_g.data(), electron_potential_.data());
    double electron_cross_energy = 0.0;
    const int density_channels = charge->nspin == 4 ? 1 : charge->nspin;
    for (int spin = 0; spin < density_channels; ++spin)
    {
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            electron_cross_energy += charge->rho[spin][ir] * electron_potential_[ir];
        }
    }
    electron_cross_energy *= dv;
#ifdef __MPI
    Parallel_Reduce::reduce_pool(electron_cross_energy);
#endif
    for (int spin = 0; spin < potential.nr; ++spin)
    {
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            potential(spin, ir) += electron_potential_[ir];
        }
    }
    energy_rydberg_ = self_energy + ion_cross_energy + electron_cross_energy;
    ModuleBase::timer::end("PotPlanarCounterCharge", "cal_v_eff");
}

void PotPlanarCounterCharge::add_force(const UnitCell&, ModuleBase::matrix& force) const
{
    const int atom_count = static_cast<int>(ion_charges_.size());
    for (int iatom = 0; iatom < atom_count; ++iatom)
    {
        const double charge = ion_charges_[iatom];
        force(iatom, 0) -= charge * ion_gradients_[iatom].x;
        force(iatom, 1) -= charge * ion_gradients_[iatom].y;
        force(iatom, 2) -= charge * ion_gradients_[iatom].z;
    }
}

} // namespace elecstate
