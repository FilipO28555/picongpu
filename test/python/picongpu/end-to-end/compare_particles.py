"""
This file is part of PIConGPU.
Copyright 2025 PIConGPU contributors
Authors: Julian Lenz
License: GPLv3+
"""

import numpy as np
import openpmd_api as opmd
import pandas as pd


def read_particles(series_name):
    series = opmd.Series(str(series_name), opmd.Access.read_only)
    names, particles = zip(*series.iterations[0].particles.items())

    data = pd.concat((particle.to_df() for particle in particles), keys=names)
    return data.assign(
        setup=data.index.to_frame()[0].apply(lambda x: x.split("_")[0]),
        impl=data.index.to_frame()[0].apply(lambda x: "_".join(x.split("_")[1:])),
    ).set_index(["setup", "impl"], drop=True)


def sort_particles(data):
    return data.groupby(["setup", "impl"]).apply(lambda df: df.sort_values(list(df.columns), axis=0))


def compare_particles_per_setup(data):
    result = True
    for setup_name, df in sort_particles(data).groupby("setup"):
        # We're doing (more than) twice the work here but that's fine:
        for lhs_name, lhs_df in df.groupby("impl"):
            for rhs_name, rhs_df in df.groupby("impl"):
                matched = np.allclose(
                    lhs_df.reset_index(drop=True),
                    rhs_df.reset_index(drop=True),
                    atol=0,
                    rtol=1.0e-4,
                )
                if not matched:
                    print(f"Mismatch found in {setup_name}: {lhs_name} != {rhs_name}")
                result *= matched
    return result


def compute_densities_per_setup_and_impl(data):
    return data.groupby(["setup", "impl", "positionOffset_x", "positionOffset_y", "positionOffset_z"])[
        "weighting"
    ].sum()


def compute_densities_from_particles(series_name):
    """
    This density is not normalised by volume yet.
    """
    return compute_densities_per_setup_and_impl(read_particles(series_name))


def _density_into_mesh(df, number_of_cells, cell_size):
    from_particles = np.zeros(number_of_cells)
    from_particles[*df.index.to_frame().to_numpy().T] = df["weighting"].to_numpy()
    return from_particles / np.prod(cell_size)


def read_densities_into_mesh(filename, number_of_cells, cell_size):
    df = (
        compute_densities_from_particles(filename)
        .reset_index(drop=False)
        .rename({"positionOffset_" + key: key for key in "xyz"}, axis=1)
    )

    for i, key in enumerate("xyz"):
        df[key] *= 1 / cell_size[i]
        df[key] = df[key].round()

    return (
        df.astype(dict(x=int, y=int, z=int))
        .set_index(["x", "y", "z"])
        .groupby(["setup", "impl"])
        .apply(
            lambda df: _density_into_mesh(df, number_of_cells, cell_size),
            include_groups=False,
        )
    )


def compare_particles(series_name):
    """
    Compare particles from the given series and return True if they all compare equal.
    """
    return compare_particles_per_setup(read_particles(series_name))


def read_binning(filename, cell_size):
    series = opmd.Series(
        str(filename),
        opmd.Access.read_only,
    )
    data = series.iterations[0].meshes["Binning"]
    mesh = data.load_chunk()
    series.flush()
    series.close()
    return mesh / np.prod(cell_size)


def read_position_check(filename):
    series = opmd.Series(
        str(filename),
        opmd.Access.read_only,
    )
    data = series.iterations[0].meshes["Binning"]
    mesh = data.load_chunk()
    series.flush()
    series.close()
    # There's an under-/overflow bin respectively.
    return mesh[1]
