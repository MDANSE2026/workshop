#!/usr/bin/env python

"""Script for processing the results of MDANSE
Dynamic Incoherent Scattering Function analysis.
It will read the results from the file, convert
the axes to the units specified in the script,
and summarise the FWHM of the central peak
as a function of Q^2.

This script assumes that MDANSE is installed in the
environment that is being used to run it. It uses
MDANSE itself as well as its dependencies."""

import numpy as np
import h5py
from scipy.interpolate import interp1d

try:
    import matplotlib.pyplot as mpl
except ImportError:
    print("Matplotlib could not be used. Plotting disabled.")
    MAKE_PLOTS = False
else:
    MAKE_PLOTS = True

from MDANSE.Framework.Units import measure


def read_data_from_mda(filename: str, dataset_name: str):
    """Extracts data from an MDA file (HDF5)

    Parameters
    ----------
    filename : str
        path to the MDA file
    dataset_name : str
        name of the dataset in the file

    Returns
    -------
    tuple(np.ndarray, dict, dict)
        a numpy array of the data, a dictionary of arrays of plotting axes,
        and a dictionary of axis units
    """
    source = h5py.File(filename)
    dataset = source[dataset_name]
    axes = dataset.attrs["axis"].split("|")
    axis_arrays, axis_units = {}, {}
    data_array = source[dataset_name][:]
    for axis_name in axes:
        axis_arrays[axis_name] = source[axis_name][:]
        axis_units[axis_name] = source[axis_name].attrs["units"]
    return data_array, axis_arrays, axis_units


def find_fwhm(xarray: np.ndarray, yarray: np.ndarray, plot=False):
    """A simple function for reading the full width at half-maximum (FWHM)
    from the plot. It literally attempts to read the width of the peak
    by interpolating the values, which makes it bad at dealing with noise
    in the results, but independent of the functional form of the peak.

    Parameters
    ----------
    xarray : np.ndarray
        the x axis of the data (e.g. energy)
    yarray : np.ndarray
        the y axis of the data (e.g. signal)
    plot : bool, optional
        if True, create a plot of the processed curve

    Returns
    -------
    float
        The FWHM of the peak
    """
    maxval = np.max(yarray)
    minval = np.min(yarray)
    halfval = (maxval + minval) / 2.0
    max_index = np.argmax(yarray)
    leftside = interp1d(yarray[: max_index + 1], xarray[: max_index + 1])
    rightside = interp1d(yarray[max_index:], xarray[max_index:])
    left_x = leftside(halfval)
    right_x = rightside(halfval)
    fwhm = np.abs(right_x - left_x)
    if plot and MAKE_PLOTS:
        selection = np.where(yarray >= 0.1 * maxval)
        mpl.plot(xarray[selection], yarray[selection], label="MDANSE S(q,$\\omega$)")
        mpl.plot([left_x, right_x], [halfval, halfval], label="FWHM")
    return fwhm


def write_results_to_file(
    filename: str, data: np.ndarray, energy_unit: str, q_unit: str
):
    """Output the DISF FWHM results into a text file.

    Parameters
    ----------
    filename : str
        name of the target file. The file will be overwritten if the file exists.
    data : np.ndarray
        an array where Q is the first column, and FWHM the second
    energy_unit : str
        unit of energy, for the header line
    q_unit : str
        unit of Q, for the header line
    """
    with open(filename, "w") as target:
        target.write(f"# Q^2[({q_unit})^2] fwhm[{energy_unit}]\n")
        for line in data:
            target.write(f"{line[0]} {line[1]}\n")


def make_fwhm_vs_qsquared(
    filename: str,
    energy_unit: str = "meV",
    q_unit: str = "1/ang",
    output_name: str = "results",
    create_partial_plots: bool = False,
):
    """The main function of this script. Takes in a name
    of the MDA file containing the DISF results,
    and produces plots and text output of the
    energy peak width as a function of the
    scattering vector length (squared).

    Parameters
    ----------
    filename : str
        Name of the MDA file (input)
    energy_unit : str, optional
        Unit to be used for energy values, by default "meV"
    q_unit : str, optional
        Unit to be used for Q vector length, by default "1/ang"
    output_name : str, optional
        Prefix used for naming files, by default "results"
    create_partial_plots : bool, optional
        If True a plot of the peak is created for each Q value, by default False
    """
    data, axes, axes_units = read_data_from_mda(filename, "/disf/s(q,f)/total")
    print(f"Data shape: {data.shape}")
    for key, value in axes.items():
        print(f"Axis {key} shape: {value.shape}")

    energy_conv_factor = measure(1.0, iunit=axes_units["disf/axes/omega"], equivalent=True).toval(
        energy_unit
    )
    q_conv_factor = measure(1.0, iunit=axes_units["disf/axes/q"], equivalent=True).toval(q_unit)

    fwhm_vs_qsq = np.zeros((data.shape[0], 2))
    for q_num, q_val in enumerate(axes["disf/axes/q"]):
        sqf = data[q_num]
        fwhm = find_fwhm(axes["disf/axes/omega"] * energy_conv_factor, sqf, plot=create_partial_plots)
        if create_partial_plots and MAKE_PLOTS:
            mpl.legend(loc=0)
            mpl.title(f"Peak width at Q={round(q_val*q_conv_factor,4)}")
            mpl.savefig(output_name + f"_peak_at_row{q_num}.png")
            mpl.close()
        fwhm_vs_qsq[q_num, 0] = (q_val * q_conv_factor) ** 2
        fwhm_vs_qsq[q_num, 1] = fwhm

    write_results_to_file(
        output_name + "_FWHM_vs_Qsq.txt", fwhm_vs_qsq, energy_unit, q_unit
    )

    if MAKE_PLOTS:
        mpl.plot(fwhm_vs_qsq[:, 0], fwhm_vs_qsq[:, 1])
        mpl.xlabel("Q${^2}$ ((" + q_unit + ")$^{2}$)")
        mpl.ylabel("S(q,$\\omega$) FWHM (" + energy_unit + ")")
        mpl.title("MDANSE results")
        mpl.savefig(output_name + "_main_plot_FWHM_vs_Qsq.png")


if __name__ == "__main__":
    # NOTE: the paths to files have to be adjusted to match the file location on your computer
    make_fwhm_vs_qsquared(
        "results/argon-disf.mda",
        energy_unit="meV",
        q_unit="1/ang",
        output_name="argon_sqw",
        create_partial_plots=True,
    )
