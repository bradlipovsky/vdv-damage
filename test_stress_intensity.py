"""Independent quadrature and physical limits of the VdV state calculation."""

import numpy as np
import pytest
from scipy.integrate import quad

from stress_intensity import stress_intensity_factors


def example_state(**changes):
    state = dict(thickness=400.0, surface_depth=20.0, basal_depth=60.0,
                 membrane_stress=[[150e3, 30e3], [30e3, 70e3]],
                 crack_normal=[np.cos(.3), np.sin(.3)], surface_water_depth=6.0,
                 basal_water_pressure=.98*917*9.81*400)
    return state | changes


@pytest.mark.parametrize("a,b,fill,support", [(4, 12, 0, 1), (80, 160, .5, 1),
                                            (280, 80, 1, .6), (20, 320, .1, .6)])
def test_against_direct_adaptive_quadrature(a, b, fill, support):
    state = example_state(surface_depth=a, basal_depth=b,
                          surface_water_depth=fill*a,
                          basal_water_pressure=support*917*9.81*400)
    h, p = state["thickness"], state["basal_water_pressure"]
    nx, ny = state["crack_normal"]
    r = 150e3*nx*nx+2*30e3*nx*ny+70e3*ny*ny

    def weight(d, z):
        return (3.52*(1-z)/(1-d)**1.5-(4.35-5.28*z)/np.sqrt(1-d)
                + ((1.30-.30*z**1.5)/np.sqrt(1-z*z)+.83-1.76*z)*(1-(1-z)*d))

    def integrate(d, traction, waterline):
        return quad(lambda z: traction(z)*weight(d, z), 0, 1,
                    points=[waterline] if 0 < waterline < 1 else None,
                    epsabs=1e-3, epsrel=1e-10)[0]

    d = a/h
    f = 1.12-.23*d+10.55*d**2-21.72*d**3+30.39*d**4
    dry = 1-fill
    ks = f*r*np.sqrt(np.pi*a)+2*np.sqrt(a/np.pi)*integrate(
        d, lambda z: -917*9.81*a*z+1000*9.81*a*max(z-dry, 0), dry)
    kb = 2*np.sqrt(b/np.pi)*integrate(
        b/h, lambda z: r-917*9.81*(h-b*z)+max(p-1024*9.81*b*z, 0), p/(1024*9.81*b))
    np.testing.assert_allclose(stress_intensity_factors(state), [ks, kb], rtol=2e-6, atol=.1)


def test_hydrostatic_cancellation_for_equal_densities():
    # A verification limit: water with the same density as ice cancels
    # hydrostatic closure in a fully filled surface crack.
    state = example_state(crack_normal=[1, 0], surface_water_depth=20)
    ks, _ = stress_intensity_factors(state, surface_water_density=917)
    d = state["surface_depth"]/state["thickness"]
    f = 1.12-.23*d+10.55*d*d-21.72*d**3+30.39*d**4
    np.testing.assert_allclose(ks, f*150e3*np.sqrt(np.pi*20), rtol=1e-12)
    # At zero resistive stress, equal-density water supporting the whole
    # column cancels the basal closure at every positive crack depth.
    state = example_state(membrane_stress=np.zeros((2, 2)), basal_depth=np.linspace(1, 350, 20),
                          basal_water_pressure=917*9.81*400)
    _, kb = stress_intensity_factors(state, basal_water_density=917)
    np.testing.assert_allclose(kb, 0, atol=1e-6)


def test_coordinate_rotation_does_not_change_intensities():
    state = example_state()
    c, s = np.cos(.7), np.sin(.7)
    q = np.array([[c, -s], [s, c]])
    rotated = state | {"membrane_stress": q@state["membrane_stress"]@q.T,
                       "crack_normal": q@state["crack_normal"]}
    np.testing.assert_allclose(stress_intensity_factors(rotated),
                               stress_intensity_factors(state), rtol=1e-12)
    reversed_normal = state | {"crack_normal": -np.array(state["crack_normal"])}
    np.testing.assert_allclose(stress_intensity_factors(reversed_normal),
                               stress_intensity_factors(state))


def test_broadcast_grid_matches_individual_columns():
    angles = np.linspace(0, np.pi, 4)
    normals = np.stack((np.cos(angles), np.sin(angles)), axis=-1)
    state = example_state(thickness=np.array([300, 400, 500])[:, None],
                          membrane_stress=np.array(example_state()["membrane_stress"])[None, None, :, :],
                          surface_depth=np.array([10, 20, 30, 40]), crack_normal=normals)
    ks, kb = stress_intensity_factors(state)
    assert ks.shape == kb.shape == (3, 4)
    for i in range(3):
        for j in range(4):
            column = state | {"thickness": state["thickness"][i, 0],
                              "surface_depth": state["surface_depth"][j],
                              "membrane_stress": state["membrane_stress"][0, 0],
                              "crack_normal": normals[j]}
            np.testing.assert_allclose([ks[i, j], kb[i, j]], stress_intensity_factors(column))


def test_water_opens_cracks_and_negative_intensity_is_preserved():
    dry = example_state(membrane_stress=np.zeros((2, 2)),
                        surface_water_depth=0, basal_water_pressure=0)
    ks, kb = stress_intensity_factors(dry)
    assert ks < 0 and kb < 0
    wet = dry | {"surface_water_depth": dry["surface_depth"],
                 "basal_water_pressure": 917*9.81*dry["thickness"]}
    ks_wet, kb_wet = stress_intensity_factors(wet)
    assert ks_wet > ks and kb_wet > kb


@pytest.mark.parametrize("changes", [dict(surface_depth=0), dict(basal_depth=380),
                                      dict(surface_water_depth=21), dict(basal_water_pressure=-1),
                                      dict(crack_normal=[2, 0]), dict(thickness=np.nan)])
def test_invalid_geometry_is_rejected_without_clipping(changes):
    with pytest.raises(ValueError):
        stress_intensity_factors(example_state(**changes))
