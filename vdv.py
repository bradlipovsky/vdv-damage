"""Van der Veen edge-crack stress intensities and Freund kinetics, in SI units.

Equations and assumptions are given in damage_law.tex. Crack fractions are
measured relative to affine vertical deformation; there is no nucleation or
healing. Full-thickness connection terminates the calculation.
"""

from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

YEAR = 365.25 * 24 * 3600
RHO_I, RHO_M, RHO_W, G = 917.0, 1000.0, 1024.0, 9.81


def slab_factor(fraction):
    """Finite-slab tension factor F(a/H), dimensionless."""
    d = np.asarray(fraction)
    return 1.12 - 0.23*d + 10.55*d**2 - 21.72*d**3 + 30.39*d**4


def slab_weight(fraction, gamma):
    """Finite-slab weight G(a/H, y/a); 0 <= gamma < 1."""
    d, z = fraction, gamma
    return (3.52*(1-z)/(1-d)**1.5 - (4.35-5.28*z)/np.sqrt(1-d)
            + ((1.30-0.30*z**1.5)/np.sqrt(1-z*z) + 0.83-1.76*z)
            * (1-(1-z)*d))


@lru_cache(maxsize=None)
def quadrature(order):
    return leggauss(order)


def weight_integral(fraction, pressure, lower=0.0, upper=1.0, order=32):
    """Integral of pressure(gamma) G(fraction,gamma) dgamma.

    gamma = sin(theta) removes the integrable crack-tip singularity. Bounds
    allow integration exactly up to a waterline, without integrating a kink.
    """
    nodes, weights = quadrature(order)
    lo, hi = np.arcsin(lower), np.arcsin(upper)
    half = np.asarray((hi-lo)/2)[..., None]
    theta = np.asarray(lo)[..., None] + half*(nodes+1)
    gamma = np.sin(theta)
    d, c = np.asarray(fraction)[..., None], np.cos(theta)
    # Algebraically cancel sqrt(1-gamma**2) against cos(theta), including
    # the zero-length interval of a dry surface crack.
    kernel = (3.52*(1-gamma)/(1-d)**1.5 - (4.35-5.28*gamma)/np.sqrt(1-d))*c
    kernel += (1.30-0.30*gamma**1.5 + (0.83-1.76*gamma)*c)*(1-(1-gamma)*d)
    return np.sum(pressure(gamma)*kernel*half*weights, axis=-1)


def stress_intensities(ds, db, thickness, stress, basal_pressure,
                       water_depth=0.0, order=32):
    """Return (K_surface, K_basal) in Pa sqrt(m).

    ds, db: positive crack fractions below 1; thickness, water_depth: m;
    stress: resistive stress normal to the cracks, Pa; basal_pressure: Pa.
    Surface water_depth is the water-column height measured from the tip.
    Both families are isolated edge cracks, as in the accompanying note.
    """
    ds, db, h, r, pb, hw = np.broadcast_arrays(
        ds, db, thickness, stress, basal_pressure, water_depth)
    if (np.any(h <= 0) or np.any(ds <= 0) or np.any(db <= 0)
            or np.any(ds >= 1) or np.any(db >= 1)):
        raise ValueError("Require H > 0 and finite flaws 0 < D_s, D_b < 1")
    a, b = h*ds, h*db
    if np.any(hw < 0) or np.any(hw > a) or np.any(pb < 0):
        raise ValueError("Require 0 <= surface water depth <= crack depth and p_b >= 0")
    ks = slab_factor(ds)*r*np.sqrt(np.pi*a)
    ks -= 2*np.sqrt(a/np.pi)*weight_integral(
        ds, lambda z: RHO_I*G*a[..., None]*z, order=order)
    dry = 1-hw/a
    ks += 2*np.sqrt(a/np.pi)*weight_integral(
        ds, lambda z: RHO_M*G*a[..., None]*(z-dry[..., None]),
        lower=dry, order=order)
    kb = 2*np.sqrt(b/np.pi)*weight_integral(
        db, lambda z: r[..., None]-RHO_I*G*(h[..., None]-b[..., None]*z),
        order=order)
    wet = np.minimum(pb/(RHO_W*G*b), 1.0)
    kb += 2*np.sqrt(b/np.pi)*weight_integral(
        db, lambda z: pb[..., None]-RHO_W*G*b[..., None]*z,
        upper=wet, order=order)
    return ks, kb


def freund_factor(q):
    """Zero below toughness, 1 - q**(-2) above; negative K cannot grow."""
    q = np.asarray(q, dtype=float)
    inverse = np.zeros_like(q)
    np.divide(1.0, q, out=inverse, where=q >= 1)
    return np.where(q >= 1, 1-inverse**2, 0.0)


def advance_damage(ds, db, thickness, stress, basal_pressure, duration,
                   speed_s, speed_b, toughness=1e5, water_fraction=0.0,
                   rtol=1e-6, atol=1e-9, order=32):
    """Integrate the local source for duration seconds, at fixed forcing.

    Speeds are m/s; toughness is Pa sqrt(m). The prescribed surface water
    column is water_fraction times the instantaneous surface-crack depth.
    Returns ds, db, elapsed seconds, and whether the first D_s+D_b=1 event
    occurred. Adaptive RK23 resolves source dynamics within a flow step.
    This array implementation is intended for a serial toy simulation.
    """
    if min(speed_s, speed_b, duration) < 0 or toughness <= 0:
        raise ValueError("Require nonnegative speeds/duration and positive toughness")
    if not 0 <= water_fraction <= 1:
        raise ValueError("Surface water fraction must lie in [0, 1]")
    ds, db, h, r, pb = np.broadcast_arrays(ds, db, thickness, stress, basal_pressure)
    shape, count = ds.shape, ds.size
    initial = np.stack((ds, db)).reshape(2, count)
    if np.any(ds+db >= 1):
        raise ValueError("The initial column has already failed")
    # Validate inputs even for arrested cracks or a zero-duration call.
    stress_intensities(ds, db, h, r, pb, water_fraction*h*ds, order)
    if duration == 0:
        return ds.copy(), db.copy(), 0.0, False
    h, r, pb = h.ravel(), r.ravel(), pb.ravel()

    def rhs(time, values):
        s, b = values.reshape(2, count)
        rate = np.zeros((2, count))
        intact = s+b < 1
        # Zero is only an extension for trial stages beyond the terminal
        # event. No physical solution is continued beyond connection.
        ks, kb = stress_intensities(
            s[intact], b[intact], h[intact], r[intact], pb[intact],
            water_fraction*h[intact]*s[intact], order)
        rate[0, intact] = speed_s/h[intact]*freund_factor(ks/toughness)
        rate[1, intact] = speed_b/h[intact]*freund_factor(kb/toughness)
        return rate.ravel()

    def failure(time, values):
        return 1-np.max(values.reshape(2, count).sum(axis=0))

    failure.terminal = True
    failure.direction = -1
    solution = solve_ivp(rhs, (0, duration), initial.ravel(), method="RK23",
                         rtol=rtol, atol=atol, events=failure)
    if not solution.success:
        raise RuntimeError(solution.message)
    s, b = solution.y[:, -1].reshape(2, *shape)
    return s, b, solution.t[-1], bool(solution.t_events[0].size)
