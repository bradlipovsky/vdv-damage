"""Surface and basal mode-I stress intensities from an ice-column state.

This is the finite-slab LEFM calculation in damage_law.tex. All inputs use
SI units. It calculates stress intensity, not crack speed or nucleation.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss


def stress_intensity_factors(state, *, ice_density=917.0,
                             surface_water_density=1000.0,
                             basal_water_density=1024.0, gravity=9.81,
                             quadrature_order=32):
    """Return (K_surface, K_basal) in Pa sqrt(m), with NumPy broadcasting.

    Required entries of the plain state dictionary:
      thickness            Ice thickness H, m.
      membrane_stress      Symmetric horizontal membrane stress, Pa;
                           trailing axes (2, 2), NOT a thickness integral.
      crack_normal         Horizontal unit normal to both crack families;
                           trailing axis (2,).
      surface_depth        Existing surface-crack penetration a_s, m.
      basal_depth          Existing basal-crack penetration a_b, m.
      surface_water_depth  Water-column height measured up from the surface
                           crack tip, m; zero is dry, a_s is fully filled.
      basal_water_pressure Water pressure at the basal crack mouth, Pa.

    M = tau_horizontal + tr(tau_horizontal) I, so n.M.n = 2*tau_nn+tau_tt.
    M must be supplied consistently with the host model's rheology; this
    function adds no damage weakening or ligament-stress amplification.
    Leading tensor/vector axes and scalar fields broadcast to the returned
    shape. Scalars produce NumPy scalars; arrays can represent a whole grid.

    Assumptions: depth-uniform resistive stress, hydrostatic ice and water,
    isolated edge cracks in a rectangular slab, and positive, finite flaws
    with a_s+a_b < H. Basal water pressure is prescribed, not inferred from
    thickness. Grounded contact and crack interactions are not modeled.
    Negative K is returned as a formal closing tendency; crack-face contact
    is not solved. Toughness belongs in a separate propagation criterion.
    """
    stress = np.asarray(state["membrane_stress"], dtype=float)
    normal = np.asarray(state["crack_normal"], dtype=float)
    if stress.shape[-2:] != (2, 2) or normal.shape[-1:] != (2,):
        raise ValueError("membrane_stress and crack_normal need trailing axes (2,2) and (2,)")
    if not np.allclose(stress, stress.swapaxes(-1, -2), rtol=1e-10, atol=1e-10):
        raise ValueError("membrane_stress must be symmetric")
    if not np.allclose(np.sum(normal**2, axis=-1), 1, rtol=1e-10, atol=1e-10):
        raise ValueError("crack_normal must be a unit vector")
    resistive_stress = np.einsum("...i,...ij,...j->...", normal, stress, normal)
    names = ("thickness", "surface_depth", "basal_depth",
             "surface_water_depth", "basal_water_pressure")
    h, a, b, water, pressure, r = np.broadcast_arrays(
        *(np.asarray(state[name], dtype=float) for name in names), resistive_stress)
    if not all(np.all(np.isfinite(value)) for value in (stress, normal, h, a, b, water, pressure)):
        raise ValueError("The ice state must contain only finite values")
    if np.any(h <= 0) or np.any(a <= 0) or np.any(b <= 0) or np.any(a+b >= h):
        raise ValueError("Require H > 0, finite flaws a_s,a_b > 0, and a_s+a_b < H")
    if np.any(water < 0) or np.any(water > a) or np.any(pressure < 0):
        raise ValueError("Require 0 <= surface_water_depth <= a_s and basal_water_pressure >= 0")
    material = np.array([ice_density, surface_water_density, basal_water_density, gravity])
    if not np.all(np.isfinite(material)) or np.any(material <= 0):
        raise ValueError("Densities and gravity must be finite and positive")
    if not isinstance(quadrature_order, (int, np.integer)) or quadrature_order < 1:
        raise ValueError("quadrature_order must be a positive integer")

    nodes, weights = leggauss(quadrature_order)

    def integral(fraction, traction, lower=0.0, upper=1.0):
        # gamma = y/a = sin(theta) removes the crack-tip singularity.
        # Multiplying G by cos(theta) analytically also handles the
        # zero-length water interval of a completely dry crack.
        lo, hi = np.arcsin(lower), np.arcsin(upper)
        half = np.asarray((hi-lo)/2)[..., None]
        theta = np.asarray(lo)[..., None]+half*(nodes+1)
        z, c, d = np.sin(theta), np.cos(theta), fraction[..., None]
        kernel = (3.52*(1-z)/(1-d)**1.5 - (4.35-5.28*z)/np.sqrt(1-d))*c
        kernel += (1.30-0.30*z**1.5+(0.83-1.76*z)*c)*(1-(1-z)*d)
        return np.sum(traction(z)*kernel*half*weights, axis=-1)

    ds, db = a/h, b/h
    factor = 1.12-0.23*ds+10.55*ds**2-21.72*ds**3+30.39*ds**4
    ks = factor*r*np.sqrt(np.pi*a)
    ks -= 2*np.sqrt(a/np.pi)*integral(ds, lambda z: ice_density*gravity*a[..., None]*z)
    dry = 1-water/a
    ks += 2*np.sqrt(a/np.pi)*integral(
        ds, lambda z: surface_water_density*gravity*a[..., None]*(z-dry[..., None]),
        lower=dry)

    kb = 2*np.sqrt(b/np.pi)*integral(
        db, lambda z: r[..., None]-ice_density*gravity*(h[..., None]-b[..., None]*z))
    # Integrate only where water pressure is positive; the pressure above
    # this waterline is zero. This is a hydrostatic boundary condition.
    wet = np.minimum(pressure/(basal_water_density*gravity*b), 1)
    kb += 2*np.sqrt(b/np.pi)*integral(
        db, lambda z: pressure[..., None]-basal_water_density*gravity*b[..., None]*z,
        upper=wet)
    return ks, kb
