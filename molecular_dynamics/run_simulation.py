#!/usr/bin/env python3
"""Run an A6B6 polymer simulation with HOOMD-blue."""

from pathlib import Path
import argparse

import gsd.hoomd
import hoomd
import hoomd.md
import numpy as np


def build_initial_frame(
    n_polymers: int,
    stickers_per_type: int,
    box_length: float,
    chain_spacing: float,
    bond_spacing: float,
    particle_mass: float,
) -> gsd.hoomd.Frame:
    particles_per_chain = 2 * stickers_per_type
    frame = gsd.hoomd.Frame()
    frame.particles.N = n_polymers * particles_per_chain
    frame.particles.types = ["A", "B"]
    frame.particles.mass = np.full(frame.particles.N, particle_mass)
    frame.configuration.box = [box_length, box_length, box_length, 0, 0, 0]

    grid_size = int(np.ceil(n_polymers ** (1 / 3)))
    positions = []
    type_ids = []
    bonds = []

    for polymer in range(n_polymers):
        ix = polymer % grid_size
        iy = (polymer // grid_size) % grid_size
        iz = polymer // (grid_size * grid_size)
        position = np.array(
            [
                (ix - (grid_size - 1) / 2) * chain_spacing,
                (iy - (grid_size - 1) / 2) * chain_spacing,
                (iz - (grid_size - 1) / 2) * chain_spacing,
            ],
            dtype=np.float64,
        )

        for particle in range(particles_per_chain):
            particle_id = polymer * particles_per_chain + particle
            positions.append(position.copy())
            type_ids.append(int(particle >= stickers_per_type))
            if particle:
                bonds.append([particle_id - 1, particle_id])
            position[0] += bond_spacing

    frame.particles.position = np.asarray(positions, dtype=np.float32)
    frame.particles.typeid = np.asarray(type_ids, dtype=np.int32)
    frame.bonds.types = ["polymer_backbone"]
    frame.bonds.group = np.asarray(bonds, dtype=np.int32)
    return frame


def make_simulation(
    initial_state: Path,
    dt: float,
    kT: float,
    gamma: float,
    attraction_strength: float,
) -> hoomd.Simulation:
    simulation = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    simulation.create_state_from_gsd(filename=str(initial_state))

    fene = hoomd.md.bond.FENEWCA()
    fene.params["polymer_backbone"] = {
        "k": 0.15,
        "r0": 10.0,
        "epsilon": 0.0,
        "sigma": 2.0,
        "delta": 0.0,
    }

    lj_nlist = hoomd.md.nlist.Cell(buffer=0.4)
    lj_nlist.exclusions = ["bond"]
    lj = hoomd.md.pair.LJ(nlist=lj_nlist, mode="shift")
    for pair in [("A", "A"), ("B", "B")]:
        lj.params[pair] = {"epsilon": 1.0, "sigma": 2.0}
        lj.r_cut[pair] = 2 ** (1 / 6) * 2.0
    lj.params[("A", "B")] = {"epsilon": 0.0, "sigma": 2.0}
    lj.r_cut[("A", "B")] = 0.0

    table_nlist = hoomd.md.nlist.Cell(buffer=0.4)
    table_nlist.exclusions = ["bond"]
    table = hoomd.md.pair.Table(nlist=table_nlist)

    r_min = 0.001
    r_cut = 1.0
    r_values = np.linspace(r_min, r_cut, 1000, endpoint=False)
    potential = -0.5 * attraction_strength * (
        1.0 + np.cos(np.pi * r_values / r_cut)
    )
    force = -0.5 * attraction_strength * (
        np.pi / r_cut
    ) * np.sin(np.pi * r_values / r_cut)
    table.params[("A", "B")] = {
        "r_min": r_min,
        "U": potential,
        "F": force,
    }
    table.r_cut[("A", "B")] = r_cut

    # Explicitly disable the same-type table interactions. A zero cutoff
    # requires a matching zero-distance, one-point table.
    for pair in [("A", "A"), ("B", "B")]:
        table.params[pair] = {
            "r_min": 0.0,
            "U": np.zeros(1),
            "F": np.zeros(1),
        }
        table.r_cut[pair] = 0.0

    integrator = hoomd.md.Integrator(dt=dt)
    integrator.forces.extend([fene, lj, table])
    langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT)
    langevin.gamma.default = gamma
    integrator.methods.append(langevin)
    simulation.operations.integrator = integrator
    return simulation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=250_000)
    parser.add_argument("--output-dir", type=Path, default=Path("md_output"))
    parser.add_argument("--write-period", type=int, default=250)
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    initial_state = output_dir / "init.gsd"
    trajectory = output_dir / "trajectory.gsd"

    frame = build_initial_frame(
        n_polymers=1000,
        stickers_per_type=6,
        box_length=286.0,
        chain_spacing=25.0,
        bond_spacing=1.0,
        particle_mass=1.0,
    )
    with gsd.hoomd.open(name=str(initial_state), mode="w") as writer:
        writer.append(frame)

    simulation = make_simulation(
        initial_state=initial_state,
        dt=0.0001,
        kT=1.0,
        gamma=0.1,
        attraction_strength=14.0,
    )
    gsd_writer = hoomd.write.GSD(
        filename=str(trajectory),
        trigger=hoomd.trigger.Periodic(args.write_period),
        mode="wb",
    )
    simulation.operations.writers.append(gsd_writer)

    print(f"Running {args.steps} steps...")
    simulation.run(args.steps)
    print(f"Wrote {trajectory}")


if __name__ == "__main__":
    main()
