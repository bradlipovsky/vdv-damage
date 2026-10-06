"""Checks against independent quadrature, LEFM arrest, and a terminal event."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from vdv import (G, RHO_I, RHO_M, RHO_W, YEAR, advance_damage,
                 freund_factor, slab_factor, slab_weight, stress_intensities)


@pytest.mark.parametrize("ds,db,water", [(0.01, 0.03, 0), (0.2, 0.4, 0.5),
                                       (0.7, 0.2, 1), (0.05, 0.8, 0)])
def test_intensities_against_adaptive_quadrature(ds, db, water):
    h, r, pb = 400.0, 130e3, 0.6*RHO_I*G*400
    a, b, delta = h*ds, h*db, h*ds*(1-water)
    # Direct integration in y, independent of the transformed Gauss rule.
    integral = lambda f, depth, d: quad(
        lambda y: f(y)*slab_weight(d, y/depth), 0, depth,
        epsabs=0.01, epsrel=1e-10, points=[delta] if depth == a else None)[0]
    ks = slab_factor(ds)*r*np.sqrt(np.pi*a)
    ks += 2/np.sqrt(np.pi*a)*integral(
        lambda y: -RHO_I*G*y + RHO_M*G*max(y-delta, 0), a, ds)
    kb = 2/np.sqrt(np.pi*b)*integral(
        lambda y: r-RHO_I*G*(h-y)+max(pb-RHO_W*G*y, 0), b, db)
    calculated = stress_intensities(ds, db, h, r, pb, water*a)
    np.testing.assert_allclose(calculated, [ks, kb], rtol=2e-6, atol=0.1)


def test_freund_threshold_and_compression():
    np.testing.assert_allclose(freund_factor([-2, 0, 0.9, 1, 2, 1e8]),
                               [0, 0, 0, 0, 0.75, 1])


def test_arrest_at_lefm_root_without_healing():
    h, r, toughness = 400.0, 150e3, 1e5
    criterion = lambda d: stress_intensities(d, .01, h, r, 0)[0]-toughness
    root = brentq(criterion, .01, .2)
    s, b, elapsed, failed = advance_damage(
        .01, .01, h, r, 0, 20*YEAR, 20/YEAR, 0, toughness, rtol=1e-8, atol=1e-11)
    assert not failed and elapsed == 20*YEAR
    assert abs(s-root) < 2e-6
    assert b == .01
    # Cracks deeper than the stable arrest point do not shrink.
    s, b, _, _ = advance_damage(root+.01, .01, h, r, 0, YEAR, 20/YEAR, 0)
    np.testing.assert_allclose([s, b], [root+.01, .01], atol=1e-12)


def test_both_components_stop_at_first_connection():
    args = (.05, .05, 100, 1e6, RHO_I*G*100, 200, 1, 1)
    coarse = advance_damage(*args, water_fraction=1, rtol=1e-6)
    fine = advance_damage(*args, water_fraction=1, rtol=1e-8, atol=1e-11)
    for s, b, time, failed in (coarse, fine):
        assert failed and 0 < time < 200
        assert s > .05 and b > .05
        assert abs(s+b-1) < 1e-10
    assert abs(coarse[2]-fine[2])/fine[2] < 1e-4


def test_requires_existing_flaws_and_an_intact_column():
    with pytest.raises(ValueError, match="finite flaws"):
        stress_intensities(0, .01, 400, 1e5, 0)
    with pytest.raises(ValueError, match="already failed"):
        advance_damage(.6, .4, 400, 1e5, 0, YEAR, 1, 1)
