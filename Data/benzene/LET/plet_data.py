import h5py
import scipp as sc

q = {
    1.24: sc.arange(
        "q", 0.03175, 1.70968 + (0.00675 / 2), 0.00675, unit=sc.Unit("1/angstrom")
    ),
    1.97: sc.arange(
        "q", 0.04001, 2.15495 + (0.00851 / 2), 0.00851, unit=sc.Unit("1/angstrom")
    ),
    3.60: sc.arange(
        "q", 0.05409, 2.91309 + (0.0115 / 2), 0.0115, unit=sc.Unit("1/angstrom")
    ),
    8.61: sc.arange(
        "q", 0.08365, 4.50511 + (0.01779 / 2), 0.01779, unit=sc.Unit("1/angstrom")
    ),
    3.62: sc.arange(
        "q", 0.13002, 2.88972 + (0.0115 / 2), 0.0115, unit=sc.Unit("1/angstrom")
    ),  # only vanadium
}

omega = {
    1.24: sc.arange("omega", -0.992, 0.992 + (0.0062 / 2), 0.0062, unit=sc.Unit("meV")),
    1.97: sc.arange(
        "omega", -1.576, 1.576 + (0.00985 / 2), 0.00985, unit=sc.Unit("meV")
    ),
    3.60: sc.arange("omega", -2.88, 2.88 + (0.018 / 2), 0.018, unit=sc.Unit("meV")),
    8.61: sc.arange("omega", -6.88, 6.88 + (0.04305 / 2), 0.04305, unit=sc.Unit("meV")),
    3.62: sc.arange(
        "omega", -2.7, 2.7 + (0.018 / 2), 0.018, unit=sc.Unit("meV")
    ),  # only vanadium
}


def load_plet_data(filename, incident_energy):
    chosen_q = q[incident_energy]
    chosen_omega = omega[incident_energy]

    file = h5py.File(filename, "r")
    id = list(file.keys())[0]

    data = sc.DataArray(
        data=sc.array(
            dims=["q", "omega"],
            values=file[id]["data"]["data"][:],
            variances=file[id]["data"]["error"][:] ** 2,
        ),
        coords={"q": chosen_q, "omega": chosen_omega},
    )
    return data
