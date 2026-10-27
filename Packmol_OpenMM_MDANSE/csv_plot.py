import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EQUIL_CSV = "equilibration.csv"
PROD_CSV = "production.csv"

# End time (ps) of the NPT, NVT, NVE-equilibration and NVE-production stages. These must
# match the stage lengths in run_water_box.py.
STAGE_ENDS = (100.0, 150.0, 175.0, 275.0)
STAGE_NAMES = ("NPT", "NVT", "NVE\nequil.", "NVE production")


def read_state_csv(path):
    """Read an OpenMM StateDataReporter CSV and clean up the column names."""
    df = pd.read_csv(path)
    df.columns = [c.replace("#", "").replace('"', "").strip() for c in df.columns]
    return df


def get(df, prefix):
    """Return the column whose name starts with prefix, or None if it is absent."""
    for c in df.columns:
        if c.lower().startswith(prefix.lower()):
            return df[c].to_numpy(dtype=float)
    return None


def shade_stages(ax, ends, label=False):
    """Shade alternate stages and optionally label them along the top of the axis."""
    starts = (0.0,) + tuple(ends[:-1])
    for i, (a, b, name) in enumerate(zip(starts, ends, STAGE_NAMES)):
        if i % 2 == 0:
            ax.axvspan(a, b, color="0.93", zorder=0)
        if label:
            ax.text((a + b) / 2, 1.02, name, transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=8)
    ax.set_xlim(0, ends[-1])


def main():
    ends = STAGE_ENDS
    production_start = ends[2]  # the production CSV restarts its clock at 0 ps

    time, quantities = [], {"Temperature": [], "Density": [], "Potential": [], "Total": []}
    for path, offset in ((EQUIL_CSV, 0.0), (PROD_CSV, production_start)):
        df = read_state_csv(path)
        t = get(df, "Time") + offset
        time.append(t)
        for key in quantities:
            quantities[key].append(get(df, key))

    t = np.concatenate(time)
    data = {k: np.concatenate(v) for k, v in quantities.items()}

    panels = [
        ("Temperature (K)", "Temperature"),
        ("Density (g/mL)", "Density"),
        ("Potential energy (kJ/mol)", "Potential"),
        ("Total energy (kJ/mol)", "Total"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    for k, (ax, (label, key)) in enumerate(zip(axes.ravel(), panels)):
        shade_stages(ax, ends, label=(k < 2))
        ax.plot(t, data[key], lw=1.0)
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("Time (ps)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()