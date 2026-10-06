"""Small, serial Icepack shelf experiment; lengths m, flow times yr, stress MPa.

The VdV routines use SI, so conversions occur explicitly at their boundary.
Run and plot the experiment with ice_shelf_damage.ipynb.
"""

import numpy as np
import firedrake as fd
import icepack
from icepack.models.viscosity import viscosity_depth_averaged, membrane_stress
from vdv import YEAR, RHO_I, G, advance_damage


class AdvectionOnly(icepack.models.DamageTransport):
    """Retain Icepack's transport; integrate the stiff VdV source separately."""

    def sources(self, **kwargs):
        return fd.Constant(0.0)


def damaged_viscosity(**fields):
    """Load-bearing-fraction closure, equivalent to A_D = A/(1-D)**n."""
    return (1-fields["damage"])*viscosity_depth_averaged(**fields)


def run_shelf(nx=40, ny=20, final_time=30.0, dt=0.25, speed=5.0,
              snapshot_interval=1.0, reaction_rtol=1e-6, progress=False):
    """Evolve a 40 x 20 km floating shelf, stopping at first connection.

    speed is the effective limiting speed of each crack family in m/yr.
    Dry surface cracks and ocean-connected basal cracks have normal (1, 0).
    A fixed front exports ice. Inflow speed is 800 sin(pi*y/W)**2 m/yr;
    lateral walls are fixed, thickness inflow is 400 m, and mass balance is
    zero. Scalar damage fractions are transported relative to affine column
    stretching. No surface/basal melting, crack rotation, or healing is used.
    """
    length, width = 40e3, 20e3
    mesh = fd.RectangleMesh(nx, ny, length, width)
    if mesh.comm.size != 1:
        raise ValueError("This toy uses serial cell arrays; run on one MPI rank")
    x, y = fd.SpatialCoordinate(mesh)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    V = fd.VectorFunctionSpace(mesh, "CG", 2)
    C = fd.FunctionSpace(mesh, "DG", 0)
    h = fd.Function(Q, name="thickness").interpolate(400-100*x/length)
    h_in = fd.Constant(400.0)
    u = fd.Function(V, name="velocity").interpolate(
        fd.as_vector((800*fd.sin(np.pi*y/width)**2, 0.0)))
    ds = fd.Function(C, name="surface damage").assign(0.01)
    patch = fd.exp(-((x-8e3)/2e3)**2-((y-10e3)/3e3)**2)
    db = fd.Function(C, name="basal damage").interpolate(0.01+0.06*patch)
    marker = fd.Function(C, name="passive initial patch").interpolate(patch)
    total = fd.Function(C, name="total damage").assign(ds+db)
    fluidity = fd.Constant(icepack.rate_factor(263.15))
    flow = icepack.solvers.FlowSolver(
        icepack.models.IceShelf(viscosity=damaged_viscosity),
        dirichlet_ids=[1, 3, 4],  # left inflow; lower and upper fixed walls
        diagnostic_solver_parameters={"snes_rtol": 1e-9})
    transport = [icepack.solvers.DamageSolver(AdvectionOnly()) for _ in range(3)]
    hc, stress = fd.Function(C), fd.Function(C)
    area = fd.assemble(fd.Constant(1)*fd.dx(domain=mesh))
    history, snapshots = [], []

    def diagnose():
        nonlocal u
        total.assign(ds+db)
        u = flow.diagnostic_solve(
            velocity=u, thickness=h, fluidity=fluidity, damage=total)
        hc.project(h)
        # M is the bulk membrane stress for the SAME damaged rheology used
        # in the momentum solve. M_xx = 2*tau_xx + tau_yy, not tau_xx alone.
        strain = fd.sym(fd.grad(u))
        stress.project((1-total)*membrane_stress(
            strain_rate=strain, fluidity=fluidity)[0, 0])

    def record(time, save=False):
        total.assign(ds+db)
        mean = lambda f: float(fd.assemble(f*fd.dx)/area)
        mass = float(fd.assemble(marker*fd.dx))
        center = float(fd.assemble(x*marker*fd.dx)/mass) if mass > 0 else np.nan
        history.append([time, mean(ds), mean(db), total.dat.data_ro.max(),
                        mean(h), center])
        if save:
            fields = (h, u, ds, db, marker, stress)
            snapshots.append(dict(zip(
                ("thickness", "velocity", "ds", "db", "marker", "stress"),
                [f.copy(deepcopy=True) for f in fields]), time=time))

    diagnose()
    record(0.0, save=True)
    time, next_snapshot, failed, max_courant = 0.0, snapshot_interval, False, 0.0
    while time < final_time-1e-12:
        step = min(dt, final_time-time, next_snapshot-time)
        # For a material fraction, the -D div(u) volume term cancels the
        # outgoing flux contribution, leaving incoming flux / cell area
        # as the diagonal loss rate. SSPRK3 inherits this Euler CFL bound.
        normal = fd.FacetNormal(mesh)
        test = fd.TestFunction(C)
        incoming = fd.assemble(
            (test("+")*fd.max_value(-fd.dot(u("+"), normal("+")), 0)
             + test("-")*fd.max_value(-fd.dot(u("-"), normal("-")), 0))*fd.dS
            + test*fd.max_value(-fd.dot(u, normal), 0)*fd.ds)
        cell_area = fd.assemble(test*fd.dx)
        loss = np.max(incoming.dat.data_ro/cell_area.dat.data_ro)
        step = min(step, 0.5/loss) if loss > 0 else step
        max_courant = max(max_courant, step*loss)

        s, b, elapsed, failed = advance_damage(
            ds.dat.data_ro, db.dat.data_ro, hc.dat.data_ro,
            stress.dat.data_ro*1e6, RHO_I*G*hc.dat.data_ro, step*YEAR,
            speed/YEAR, speed/YEAR, rtol=reaction_rtol)
        ds.dat.data[:] = s
        db.dat.data[:] = b
        if failed:
            # The source event is terminal; never solve flow at zero
            # load-bearing thickness or silently cap D below one.
            time += elapsed/YEAR
            record(time, save=True)
            break

        for solver, field, inflow in zip(transport, (ds, db, marker), (0.01, 0.01, 0.0)):
            field.assign(solver.solve(step, damage=field, velocity=u,
                                      damage_inflow=fd.Constant(inflow)))
        h = flow.prognostic_solve(step, thickness=h, velocity=u,
                                 accumulation=fd.Constant(0), thickness_inflow=h_in)
        if h.dat.data_ro.min() <= 0:
            raise RuntimeError("Nonpositive thickness; reduce the flow time step")
        time += step
        diagnose()
        save = time >= next_snapshot-1e-10 or time >= final_time-1e-10
        record(time, save=save)
        if save:
            next_snapshot += snapshot_interval
            if progress:
                print(f"t={time:5.1f} yr, max D={history[-1][3]:.3f}, "
                      f"mean H={history[-1][4]:.1f} m", flush=True)

    return {"mesh": mesh, "snapshots": snapshots, "history": np.array(history),
            "failed": failed, "max_courant": max_courant,
            "speed": speed, "dt": dt, "resolution": (nx, ny)}
