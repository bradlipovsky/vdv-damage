"""Transport convergence and a coupled shelf calculation; requires Firedrake."""

import numpy as np
import firedrake as fd
import icepack
from ice_shelf import AdvectionOnly, run_shelf


def test_advection_converges_to_translation():
    errors = []
    for nx in (16, 32):
        mesh = fd.PeriodicRectangleMesh(nx, nx//2, 4000, 2000, direction="both")
        x, y = fd.SpatialCoordinate(mesh)
        C = fd.FunctionSpace(mesh, "DG", 0)
        V = fd.VectorFunctionSpace(mesh, "CG", 1)
        u = fd.Function(V).assign(fd.Constant((1000, 200)))
        initial = .1+.05*fd.sin(2*np.pi*x/4000)*fd.cos(2*np.pi*y/2000)
        damage = fd.Function(C).project(initial)
        mass = fd.assemble(damage*fd.dx)
        solver = icepack.solvers.DamageSolver(AdvectionOnly())
        end = .5
        steps = int(np.ceil(end/(.25*(4000/nx)/1200)))
        for _ in range(steps):
            damage.assign(solver.solve(end/steps, damage=damage, velocity=u))
        exact = .1+.05*fd.sin(2*np.pi*(x-1000*end)/4000)*fd.cos(2*np.pi*(y-200*end)/2000)
        errors.append(np.sqrt(fd.assemble((damage-exact)**2*fd.dx)/(4000*2000)))
        assert abs(fd.assemble(damage*fd.dx)/mass-1) < 1e-10
        assert damage.dat.data_ro.min() >= .05-1e-12
        assert damage.dat.data_ro.max() <= .15+1e-12
    print("DG0 translation RMS errors:", errors)
    assert errors[1] < .8*errors[0]


def test_coupled_shelf_grows_and_transports_damage():
    result = run_shelf(nx=12, ny=6, final_time=1, snapshot_interval=.5)
    history = result["history"]
    assert not result["failed"]
    assert history[-1, 0] == 1
    assert history[-1, 2] > history[0, 2]
    assert history[-1, 5] > history[0, 5]
    assert result["max_courant"] <= .5+1e-12
    for frame in result["snapshots"]:
        assert frame["thickness"].dat.data_ro.min() > 0
        assert np.all(frame["ds"].dat.data_ro > 0)
        assert np.all(frame["db"].dat.data_ro > 0)
        assert np.max(frame["ds"].dat.data_ro+frame["db"].dat.data_ro) < 1


def test_fraction_stays_constant_in_divergent_flow():
    # D is a material fraction, not a conserved areal density.
    mesh = fd.RectangleMesh(8, 4, 4000, 2000)
    x, y = fd.SpatialCoordinate(mesh)
    C = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    u = fd.Function(V).interpolate(fd.as_vector((100+x/10, y/20)))
    damage = fd.Function(C).assign(.2)
    solver = icepack.solvers.DamageSolver(AdvectionOnly())
    for _ in range(4):
        damage.assign(solver.solve(.05, damage=damage, velocity=u,
                                   damage_inflow=fd.Constant(.2)))
    np.testing.assert_allclose(damage.dat.data_ro, .2, atol=1e-12)


def test_shelf_stops_before_solving_flow_at_full_damage():
    result = run_shelf(nx=8, ny=4, final_time=1, speed=1e4)
    assert result["failed"]
    assert 0 < result["history"][-1, 0] < .25
    assert abs(result["history"][-1, 3]-1) < 1e-9
    assert len(result["snapshots"]) == 2
    np.testing.assert_array_equal(
        result["snapshots"][0]["velocity"].dat.data_ro,
        result["snapshots"][-1]["velocity"].dat.data_ro)
