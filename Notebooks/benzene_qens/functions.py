import numpy as np
from scipy.integrate import trapezoid


def second_moment_analyser(energy, s_qw, s_qw_err=None):
    """
    Calculate the second moment of S(q,ω).

    Parameters
    ----------
    energy: ndarray
         Array of energy values.
    s_qw: ndarray
        2D array where each row corresponds to the scattering function at a specific energy.

    Returns
    ----------
        second_moment: ndarray
            Array of second moment values as a function of q
    """
    s_qe = np.zeros(s_qw.shape[0])
    e2s_qe = np.zeros(s_qw.shape[0])

    if s_qw_err is not None:
        err_total = np.zeros(s_qw.shape[0])

    for i, sqw in enumerate(s_qw):
        s_qe[i] = trapezoid(y=sqw, x=energy)
        e2s_qe[i] = trapezoid(y=(energy**2 * sqw), x=energy)
        if s_qw_err is not None:
            _, s_qe_err = trap_uncertainty(x=energy, y=sqw, var_y=s_qw_err[i])
            _, e2s_qe_err = trap_uncertainty(
                x=energy, y=(energy**2 * sqw), var_y=(energy**2 * s_qw_err[i])
            )
            err_total[i] = np.sqrt(
                (e2s_qe[i] / s_qe[i] ** 2) ** 2 * s_qe_err**2
                + (1 / s_qe[i]) ** 2 * e2s_qe_err**2
            )

    second_moment = e2s_qe / s_qe

    if s_qw_err is not None:
        return second_moment, err_total
    else:
        return second_moment


def limit_integral(energy, s_qw, limits, s_qw_err=None):
    """
    Numerical integration of S(q,ω) with limits

    Parameters
    ----------
    energy: ndarray
         Array of energy values.
    s_qw: ndarray
        2D array where each row corresponds to the scattering function at a specific energy.
    limits: ndarray
        2 limits i.e. [-0.02,0.02]

    Returns
    ----------
        integral: ndarray
            Integral within limits as a function of q
    """
    lims = np.where((limits[0] < energy) & (energy < limits[1]))[0]
    s_qw_trimmed = s_qw[:, lims]

    if s_qw_err is not None:
        s_qw_err_trimmed = s_qw_err[:, lims]
        int_err = np.zeros(s_qw_trimmed.shape[0])

    integral = np.zeros(s_qw_trimmed.shape[0])

    for i, sqw in enumerate(s_qw_trimmed):
        integral[i] = trapezoid(y=sqw, x=energy[lims])
        if s_qw_err is not None:
            _, int_err[i] = trap_uncertainty(
                x=energy[lims], y=sqw, var_y=s_qw_err_trimmed[i]
            )

    if s_qw_err is not None:
        return integral, int_err
    else:
        return integral


def trap_uncertainty(x, y, var_y):
    """
    Compute trapezoidal integral and propagated uncertainty for unevenly spaced data.
    Based off:
    https://stats.stackexchange.com/questions/214850/propagate-errors-in-measured-points-to-simpsons-numerical-integral

    Parameters
    ----------
    x : array-like
        The x-values (must be in increasing order).
    y : array-like
        The corresponding y-values.
    var_y : array-like
        variances of y-values.

    Returns
    -------
    integral : float
        The computed trapezoidal integral.
    sigma_integral : float
        The propagated uncertainty in the integral.
    """
    N = len(x)
    h = np.diff(x)

    # Trapezoidal integral
    integral = np.sum(0.5 * h * (y[:-1] + y[1:]))

    # Error propagation
    sigma_I_sq = 0.0
    sigma_I_sq += (0.5 * h[0]) ** 2 * var_y[0]
    sigma_I_sq += (0.5 * h[-1]) ** 2 * var_y[-1]

    for j in range(1, N - 1):
        coeff = 0.5 * (h[j - 1] + h[j])
        sigma_I_sq += coeff**2 * var_y[j]

    sigma_integral = np.sqrt(sigma_I_sq)

    # Assert solution is similar to scipy implementation
    assert np.isclose(integral, trapezoid(x=x, y=y))

    return integral, sigma_integral
