#! /usr/bin/env python3
"""
Script to compute nonlinear transfer function

Copyright 2020 Steve Biggs

This file is part of initial.

initial is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

initial is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with initial.  If not, see http://www.gnu.org/licenses/.

Modified by Tobias Schuett on 17 June 2024
"""

import multiprocessing as mp
import numpy as np
import sys
import xarray as xr


def parse_args():
    input_filename, it_start, it_end, nproc = sys.argv[1:]
    it_start = int(it_start)
    it_end = int(it_end)
    nproc = int(nproc)
    return input_filename, it_start, it_end, nproc


def parallel_compute_transfer_along_theta():
    # Loop over theta in parallel
    with mp.Pool(processes=nproc) as p:
        result = np.array(
            p.map(parallel_compute_transfer_wrapper, range(ntheta), ntheta // nproc),
            dtype=float,
        )
    return result


def parallel_compute_transfer_wrapper(i_theta):
    print("Computing net transfer for theta = {}".format(ds["theta"].values[i_theta]))
    return compute_net_transfer_this_theta(
        ds["phi_t"][:, :, :, i_theta, :], ds["bmag"].values[i_theta],
    )


def compute_net_transfer_this_theta(phi, B):
    # Convert ri dimension to complex dtype
    phi = ri_to_complex(phi)
    # Move t axis to last position for improved cache performance
    phi = move_first_axis_to_last(phi)
    # Make full versions of arrays to simplify indexing
    phi = make_full(phi)
    # Pre-shift arrays so that indexing works correctly
    phi = np.fft.fftshift(phi, axes=(0, 1))
    # Pre-compute complex conjugate array
    phi_star = np.conj(phi[ikyt, ...])
    # Pre-prepare array of mediators for performance
    phi_m = compute_phi_m(phi, ikx0, iky0)
    # calculate B_factor
    B_factor = 1 / B**3
    # Compute transfer and return
    return compute_net_transfer(
        phi, phi_star, phi_m, B_factor
    )


def ri_to_complex(array_ri):
    array_complex = np.zeros(array_ri.shape[:-1], dtype=complex)
    array_complex.real = array_ri[..., 0]
    array_complex.imag = array_ri[..., 1]
    return array_complex


def move_first_axis_to_last(array_orig):
    # Roll first axis to last place, e.g. for improved cache performance
    # NB: Cannot use np.moveaxis(...) as the return value is just a view of the original
    array_new = np.zeros(
        np.concatenate((array_orig.shape[1:], [array_orig.shape[0]])),
        dtype=array_orig.dtype,
    )
    for i in range(array_orig.shape[0]):
        array_new[..., i] = array_orig[i, ...]
    return array_new


def make_full(array):
    # Assumes dimension order is ky, kx, t
    # For concatenated array:
    # - ky slice = -1:0:-1 so that we include non-zero kys in reverse order
    # - kx slice = -1::-1 so that we include all kxs in reverse order but this puts kx = 0
    #     to the end even though we still want it at the start, hence use of
    #     np.roll(..., 1, axis=1) so that kx = 0 is at ikx = 0
    # - np.conj(...) of the reversed and rolled array as the negative kys are the conjugate
    #     of the positive kys
    # - np.concatenate(..., axis=0) so that we concatenate in the ky direction
    return np.concatenate(
        (array, np.conj(np.roll(array[-1:0:-1, -1::-1, :], 1, axis=1))),
        axis=0,
    )


def compute_phi_m(phi, ikx0, iky0):
    # Pre-allocate array
    # NB: Cannot pre-allocate this variable in the global scope and re-use it here
    phi_m = np.zeros((nky, nkx, nkx, nt), dtype=complex)
    # Loop over target and source wavenumbers
    for ikys in range(nky):
        for ikxs in range(nkx):
            for ikxt in range(nkx):
                # Work out index of mediator
                ikxm = ikxt - ikxs + ikx0
                ikym = ikyt - ikys + iky0
                # Check mediator index exists
                if not (0 <= ikxm and ikxm < nkx and 0 <= ikym and ikym < nky):
                    # Just don't set a value to avoid unnecessary cache misses
                    continue
                # Store mediator value in mediator array
                phi_m[ikys, ikxs, ikxt, :] = phi[ikym, ikxm, :]
    # Return output array
    return phi_m


def compute_net_transfer(phi, phi_star, phi_m, B_factor):
    T_v = (
        2 * B_factor
        * z_hat_dot_k_cross_k_prime
        * k_factor_for_T_v
        * compute_time_average(
            (
                np.reshape(phi_star, (1, 1, nkx, nt))
                * phi_m
                * np.reshape(phi, (nky, nkx, 1, nt))
            ).real
        )
    )

    return np.sum(T_v, axis=(0,2)) # shape (kys, kxs, kxt, t) -> (kxs, t)


def compute_time_average(field):
    "Computes the time average of an N-D field assuming time is the last axis"
    if len(t) == 1:
        print("len(t) == 1 -> not actually computing time average, just taking out time dimension")
        time_average = field[..., 0]
    else:
        print(f"len(t) == {len(t)} so actually computing time average")
        print("Computing time average...")
        time_average = (
            np.trapz(
                field,
                t,
            )
            / (t[-1] - t[0])
        )
    return time_average

def write_output():
    np.savez("Tv_theta_kxsource.npz", theta=theta, kxs=kx, T_v=result)


if __name__ == "__main__":
    # Parse command line inputs
    input_filename, it_start, it_end, nproc = parse_args()
    # Open Dataset
    if "," in input_filename:
        input_filename_list = [f for f in input_filename.split(",")]
        ds = xr.open_mfdataset(input_filename_list)
    else:
        ds = xr.open_dataset(input_filename)
    # Remove un-used variables
    required_vars = ["theta", "kx", "ky", "t", "phi_t", "bmag"]
    #vars_to_delete = ["egrid", "lambda"]  # un-used coords
    vars_to_delete = ["lambda"]  # changed for new diagnotics output

    for v in ds.data_vars:
        if v not in required_vars:
            vars_to_delete.append(v)
    ds = ds.drop_vars(vars_to_delete)

    ## NOTE: negative indexing -> positive indexing #######
    print("it_start: ", it_start)
    print("it_end: ", it_end)
    if it_start < 0:
        it_start = len(ds["t"]) + it_start
    if it_end < 0:
        it_end = len(ds["t"]) + it_end
    print("len(t): ", len(ds["t"]))
    print("it_start: ", it_start)
    print("it_end: ", it_end)
    #######################################################

    ds = ds.drop_sel(t=ds["t"][it_end:])
    ds = ds.drop_sel(t=ds["t"][:it_start])
    # Load data for performance and to avoid concurrent read problems
    ds.load()
    # Pre-compute various quantities for performance
    theta = ds["theta"].values # added TMS
    ntheta = len(theta)
    kx = np.fft.fftshift(ds["kx"].values)
    nkx = len(kx)
    ky = np.fft.fftshift(np.concatenate((ds["ky"].values, -ds["ky"].values[-1:0:-1])))
    nky = len(ky)
    ikx0 = np.argmin(np.abs(kx))
    iky0 = np.argmin(np.abs(ky))
    ikyt = iky0  # i.e. target is always zonal flows
    t = ds["t"].values
    nt = len(t)
    z_hat_dot_k_cross_k_prime = np.reshape(
        np.reshape(kx, (1, nkx)) * np.reshape(ky, (nky, 1)),
        (nky, 1, nkx),
    )
    k_factor_for_T_v = -np.reshape(
        np.reshape(kx, (1, nkx)) * np.reshape(kx, (nkx, 1)),
        (1, nkx, nkx),
    )
    # Do calculation
    result = parallel_compute_transfer_along_theta()
    # Write output
    write_output()