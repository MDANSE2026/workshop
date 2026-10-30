import h5py
import scipp as sc
import numpy as np

q_values = np.array(
    [
        0.44168133382170255,
        0.4836718818802158,
        0.5252494850177617,
        0.5667096154060947,
        0.6078712434417826,
        0.6487126881042931,
        0.689059119666155,
        0.7291972609506315,
        0.7688009107323014,
        0.8081515132390025,
        0.8470764389591491,
        0.8854097486319593,
        0.9234238587847176,
        0.960809830116343,
        0.9978333468728936,
        1.0343312759325145,
        1.0702843927798718,
        1.105540290816107,
        1.1403495168734994,
        1.1745580878433999,
        1.2080214609561863,
        1.2409774444254584,
        1.2731582061631168,
        1.3047924394761778,
        1.3357394008561196,
        1.365982789642908,
        1.395395800040148,
        1.42418745626771,
        1.4521237695746756,
        1.4794032699998334,
        1.505903525914927,
        1.5316105788683587,
        1.5564177547948364,
        1.580501346447981,
        1.6037524423088392,
        1.6260752149772537,
        1.6476283017392186,
        1.6683135349991718,
        1.6880463869882285,
        1.7069670896628295,
        1.724921875155287,
        1.742038339552496,
        1.7582372219489453,
        1.7735099899245481,
        1.7877958583340137,
        1.8011471901015434,
        1.8136025596776764,
        1.82510265230754,
        1.8356414105595806,
        1.8451787261224841,
        1.8537823761437875,
    ]
)


def rebin_density(da, **bins):
    # scipp's rebin conserves counts, but S(Q,w) is a density (per meV) on a
    # non-uniform energy axis -- convert to counts, rebin, convert back.
    out = da.copy()
    for dim, target in bins.items():
        w = sc.array(
            dims=[dim],
            values=np.diff(out.coords[dim].values),
            unit=out.coords[dim].unit,
        )
        nw = sc.array(dims=[dim], values=np.diff(target.values), unit=target.unit)
        out = (out * w).rebin({dim: target}) / nw
    return out


def load_iris(file):
    with h5py.File(file, "r") as f:
        ws = f["mantid_workspace_1/workspace"]
        I = ws["values"][...]
        dI = ws["errors"][...]
        omega = ws["axis1"][...]  # bin edges, meV
        q = q_values[ws["axis2"][...].astype(int)]

    sqw = sc.DataArray(
        data=sc.array(
            dims=["q", "omega"], values=I, variances=dI**2, unit="angstrom/meV"
        ),
        coords={
            "q": sc.array(dims=["q"], values=q, unit="1/angstrom"),
            "omega": sc.array(dims=["omega"], values=omega, unit="meV"),
        },
    )

    # dividing by the detector count per bin averages rather than sums
    q_bins = sc.linspace("q", 0.42, 1.86, num=31, unit="1/angstrom")
    ones = sc.DataArray(sc.ones(sizes=sqw.sizes), coords={"q": sqw.coords["q"]})
    binned = (sqw.hist(q=q_bins) / ones.hist(q=q_bins)).transpose(["q", "omega"])

    # uniform energy grid at the coarsest native width
    step = np.diff(omega).max()
    num = int((omega[-1] - omega[0]) // step)
    omega_bins = sc.linspace("omega", omega[0], omega[-1], num=num + 1, unit="meV")
    return rebin_density(binned, omega=omega_bins)
