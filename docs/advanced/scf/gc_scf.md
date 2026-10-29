# Grand-canonical SCF prototype

ABACUS provides a plane-wave prototype for SCF calculations at a fixed work
function. The electron number is not fixed during the SCF loop. Instead, the
Fermi-Dirac occupations are generated from a chemical potential, and that
chemical potential is mixed toward the value required by the target work
function.

## Restrictions

The current prototype is intentionally limited:

- `basis_type pw`
- `esolver_type ksdft`
- `calculation scf`
- `smearing_method fd`
- `scf_thr_type 1`
- `assume_isolated pcc_2d`
- `nspin 1` or `2`
- no fixed spin difference, DFT+U, DFT-1/2, electric fields, gate fields, or
  force/stress evaluation

The input must provide enough `nbands` for the electron-number fluctuation.

## Example

```text
INPUT_PARAMETERS
  basis_type pw
  calculation scf
  smearing_method fd
  assume_isolated pcc_2d
  pcc_2d_axis 2
  scf_thr_type 1
  nbands 36
  gce true
  gce_target_workfunction 3.0
  gce_mixing_beta 0.2
  gce_wf_tol 5.0e-3
  gce_nelec_tol 1.0e-5
  gce_qk 3.0
  gce_qkappa 0.25
  gce_qm 0.75
```

The prototype uses the planar-averaged electrostatic potential at the vacuum
density minimum as the work-function reference. It reports the current work
function, electron number, charged electron count, and the electron-reservoir
free-energy correction in the SCF log.

The PCC 2D correction supplies the charged-slab electrostatic reference. It is
not the Gaussian planar counter-charge model used by the corresponding CP2K
implementation.

`gce_qk`, `gce_qkappa`, and `gce_qm` control the modified Kerker preconditioner
and finite-`G=0` Pulay metric used by the GC density mixing.
