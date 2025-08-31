import os
from itertools import product
from pathlib import Path
import pickle

import joblib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')  # Use TkAgg backend for interactive plots

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from scipy.interpolate import griddata

def _find_path_to_data(lattice_cell):
    """Determine a default dataset path based on lattice geometry."""
    path_dataset = Path(__file__).parents[2] / "data" / "outputs" / "relative_densities" / "data"
    geom_cell = lattice_cell.cells[0].geom_types
    name_dataset = "RelativeDensities" + "".join(f"_{g}" for g in geom_cell)
    return path_dataset, name_dataset

def compute_relative_densities_dataset(lattice_cell,
                                       step_radius: float = 0.01,
                                       range_radius: tuple = (0.00, 0.1),
                                       name_dataset: str = None,
                                       save_every: int = 1,
                                       resume: bool = True):
    """
    Compute relative densities for a range of radii in the lattice.
    Periodically saves partial results every `save_every` iterations and can resume if a file exists.
    """
    if lattice_cell.get_number_cells() != 1:
        raise ValueError("The lattice must contain exactly one cell.")

    path_dataset, name_dataset = _find_path_to_data(lattice_cell)

    # Prepare grid and potential resume
    n_geom = len(lattice_cell.cells[0].geom_types)
    valid_combos = _valid_combinations(step_radius, range_radius, n_geom)

    relative_densities = {}
    dataset_file = path_dataset / f"{name_dataset}.pkl"
    if dataset_file.exists() and resume:
        relative_densities = load_dataset(path_dataset, name_dataset)
        print(f"Resuming from {len(relative_densities)} existing entries in {dataset_file}")

        # Use check_missing_entries to report and sanitize
        status = check_missing_entries(
            path_dataset=path_dataset,
            name_dataset=name_dataset,
            step_radius=step_radius,
            range_radius=range_radius,
            n_geom=n_geom,
            threshold=0.001,
        )

        # Remove invalid entries that don't match current grid/threshold
        if status.get("n_invalid", 0) > 0:
            print(f"Removing {status['n_invalid']} invalid entries from {dataset_file}")
            for combo in status["invalid_combinations"]:
                relative_densities.pop(combo, None)
            save_dataset(path_dataset, name_dataset, relative_densities)

        remaining = status["missing_combinations"]
        print(f"{len(remaining)} combinations missing (expected {status['n_expected']}).")
    else:
        remaining = [c for c in valid_combos if c not in relative_densities]

    if not remaining:
        print("Dataset already complete.")
        return

    for i, combo in enumerate(remaining, 1):
        print(f"Computing for radii: {combo}")
        if i <= 1:
            combo = list(combo)
            combo[0] += 0.001
            lattice_cell.change_beam_radius(combo)
            combo[0] -= 0.001
            combo = tuple(combo)
        else:
            lattice_cell.change_beam_radius(list(combo))
        rd = lattice_cell.generate_mesh_lattice_Gmsh(volume_computation=True,
                                                     cut_mesh_at_boundary=True,
                                                     save_STL=False,
                                                     only_relative_density=True)
        if rd is None:
            continue
        relative_densities[combo] = rd
        print(f"Relative density: {rd}")

        if i % max(1, int(save_every)) == 0:
            save_dataset(path_dataset, name_dataset, relative_densities)
            print(f"Progress saved: {len(relative_densities)}/{len(valid_combos)} entries")

    # Final save
    save_dataset(path_dataset, name_dataset, relative_densities)
    print(f"✅ Completed. Total entries: {len(relative_densities)}/{len(valid_combos)}")


def save_dataset(path_dataset, name_dataset, relative_densities_dict):
    """Save the dataset as a pickle file (atomic write)."""
    path_dataset.mkdir(parents=True, exist_ok=True)
    final_path = path_dataset / f"{name_dataset}.pkl"
    tmp_path = path_dataset / f"{name_dataset}.pkl.tmp"
    with open(tmp_path, "wb") as f:
        pickle.dump(relative_densities_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, final_path)  # atomic on POSIX/Windows
    print(f"Dataset saved as pickle at {final_path}")


def load_dataset(path_dataset, name_dataset):
    """
    Load a dataset previously saved as a pickle file.

    Parameters
    ----------
    path_dataset : Path
        Path to the directory containing the dataset.
    name_dataset : str
        Name of the dataset (without extension).

    Returns
    -------
    dict
        Dictionary of relative densities keyed by radii tuple.
    """
    file_path = path_dataset / f"{name_dataset}.pkl"
    if not file_path.exists():
        raise FileNotFoundError(f"Dataset file not found: {file_path}")

    with open(file_path, "rb") as f:
        relative_densities_dict = pickle.load(f)

    print(f"Dataset {name_dataset} loaded from {file_path}")
    return relative_densities_dict

def _valid_combinations(step_radius: float, range_radius: tuple, n_geom: int, threshold: float = 0.001):
    """Generate all valid radius combinations given the grid and a sum threshold."""
    min_radius, max_radius = range_radius
    radius_values = np.arange(min_radius, max_radius + step_radius, step_radius)
    all_combinations = list(product(radius_values, repeat=n_geom))
    return [c for c in all_combinations if sum(c) > threshold]

def check_missing_entries(path_dataset: Path,
                          name_dataset: str,
                          step_radius: float,
                          range_radius: tuple,
                          n_geom: int,
                          threshold: float = 0.001):
    """
    Check which radius combinations are missing in an existing dataset file (if any).
    Returns a dict with counts and the list of missing combinations.
    """
    expected = set(_valid_combinations(step_radius, range_radius, n_geom, threshold))
    print(len(expected), "expected combinations based on parameters.")
    try:
        data = load_dataset(path_dataset, name_dataset)
        present = set(data.keys())
    except FileNotFoundError:
        present = set()

    missing = sorted(expected - present)
    invalid = sorted(present - expected)

    if invalid:
        print(f"⚠️ Found {len(invalid)} invalid combinations in dataset:")
        for combo in invalid[:20]:
            print("   ", combo)
        if len(invalid) > 20:
            print(f"   ... and {len(invalid)-20} more")

    return {
        "n_expected": len(expected),
        "n_present": len(present),
        "n_missing": len(missing),
        "missing_combinations": missing,
        "n_invalid": len(invalid),
        "invalid_combinations": invalid,
    }


def plot_3D_iso_surface(lattice_cell):
    """
    Plot 3D iso-surfaces of volume as a function of three radii using interpolation.

    Parameters
    ----------
    lattice_cell : Lattice
        Lattice object with exactly one cell and three geometry types.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("Plotly is required for 3D plotting. Install it via 'conda install plotly'.")
    # --- Load dataset and build (radii, volumes) arrays ---
    path_dataset, name_dataset = _find_path_to_data(lattice_cell)
    dataset_file = path_dataset / f"{name_dataset}.pkl"
    if not dataset_file.exists():
        raise FileNotFoundError(f"No dataset found at {dataset_file}")

    relative_densities = load_dataset(path_dataset, name_dataset)
    if not relative_densities:
        raise ValueError("Loaded dataset is empty.")

    radii = np.array(list(relative_densities.keys()), dtype=float)    # shape (N, n_geom)
    volumes = np.array(list(relative_densities.values()), dtype=float)  # shape (N,)

    if radii.ndim != 2 or radii.shape[1] != 3:
        raise ValueError(f"Expected 3 radii per sample, got shape {radii.shape}")

    # --- Build regular 3D grid over observed radii domain ---
    grid_x, grid_y, grid_z = np.mgrid[
                             np.min(radii[:, 0]):np.max(radii[:, 0]):30j,
                             np.min(radii[:, 1]):np.max(radii[:, 1]):30j,
                             np.min(radii[:, 2]):np.max(radii[:, 2]):30j
                             ]

    # --- Interpolate volumes on the grid ---
    grid_vol = griddata(radii, volumes, (grid_x, grid_y, grid_z), method='linear')

    fig = go.Figure(data=go.Isosurface(
        x=grid_x.flatten(), y=grid_y.flatten(), z=grid_z.flatten(),
        value=grid_vol.flatten(),
        opacity=0.3, surface_count=10, colorscale='viridis'
    ))

    fig.update_layout(
        scene=dict(xaxis_title='Radius 1', yaxis_title='Radius 2', zaxis_title='Radius 3'),
        title="3D Iso-Surface of Relative Densities"
    )

    fig.show()




def evaluate_kriging_from_pickle(
    dataset_dir: Path,
    name_dataset: str,
    test_size: float = 0.2,
    model_name: str = "kriging_model_",
    kernel: object | None = None,
    normalize_y: bool = True,
    random_state: int = 42,
):
    """
    Train and evaluate a Kriging (GPR) model from a relative-density dataset saved as a pickle.

    Parameters
    ----------
    dataset_dir : Path
        Directory that contains '<name_dataset>.pkl' produced by `compute_relative_densities_dataset`.
    name_dataset : str
        Dataset base name (without extension).
    test_size : float, optional
        Fraction of samples used for testing (default 0.2).
    model_name : str | Path, optional
        Relative path (from the project root) where the trained model will be saved.
    kernel : sklearn.gaussian_process.kernels.Kernel | None, optional
        Custom GPR kernel. If None, a sensible default is used.
    normalize_y : bool, optional
        Whether to normalize the target inside the regressor (default True).
    random_state : int, optional
        Random seed for reproducibility.

    Returns
    -------
    dict
        Evaluation metrics and paths. Keys: 'MSE', 'RMSE', 'NRMSE', 'MAE', 'R2',
        'n_train', 'n_test', 'model_path', 'kernel_'
    """
    # Load dataset
    pkl_path = Path(dataset_dir) / f"{name_dataset}.pkl"
    if not pkl_path.exists():
        raise FileNotFoundError(f"Dataset file not found: {pkl_path}")

    with open(pkl_path, "rb") as f:
        rel_dens_dict = pickle.load(f)

    if not isinstance(rel_dens_dict, dict) or len(rel_dens_dict) == 0:
        raise ValueError("Loaded dataset is empty or not a dictionary.")

    # Convert dict to arrays
    try:
        X = np.array(list(rel_dens_dict.keys()), dtype=float)
        y = np.array(list(rel_dens_dict.values()), dtype=float)
    except Exception as e:
        raise ValueError(f"Failed to parse dataset content into arrays: {e}")

    if X.ndim == 1:
        X = X.reshape(-1, 1)

    # Sanity checks
    if np.any(~np.isfinite(X)) or np.any(~np.isfinite(y)):
        raise ValueError("Non-finite values detected in dataset (NaN/Inf).")

    if np.ptp(y) == 0:
        raise ValueError("All relative densities are identical; cannot evaluate NRMSE/R² meaningfully.")

    # Train / test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=float(test_size), random_state=random_state
    )

    # Kernel definition
    if kernel is None:
        kernel = (C(1.0, (1e-3, 1e3)) *
                  RBF(length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-3, 1e3))
                  + WhiteKernel(noise_level=1e-6, noise_level_bounds=(1e-12, 1e-1)))

    gpr = GaussianProcessRegressor(
        kernel=kernel,
        n_restarts_optimizer=10,
        normalize_y=normalize_y,
        random_state=random_state,
    )

    # Fit
    gpr.fit(X_train, y_train)

    # Predict & metrics
    y_pred, y_std = gpr.predict(X_test, return_std=True)

    mse = mean_squared_error(y_test, y_pred)
    rmse = float(np.sqrt(mse))
    mae = mean_absolute_error(y_test, y_pred)
    nrmse = float(rmse / (np.max(y_test) - np.min(y_test)))
    r2 = r2_score(y_test, y_pred)

    print("✅ Kriging model evaluation")
    print(f"   • Train size = {len(y_train)}, Test size = {len(y_test)}")
    print(f"   • MSE   = {mse:.6e}")
    print(f"   • RMSE  = {rmse:.6e}")
    print(f"   • MAE   = {mae:.6e}")
    print(f"   • NRMSE = {nrmse:.6e}")
    print(f"   • R²    = {r2:.6f}")
    print(f"   • Learned kernel: {gpr.kernel_}")

    # Save model
    project_root = Path(__file__).resolve().parents[2] / "data" / "outputs" / "relative_densities" / "surrogate_model"
    model_path = (project_root / (model_name + name_dataset)).resolve()
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": gpr,
            "kernel": gpr.kernel_,
            "feature_dim": X.shape[1],
            "metadata": {
                "dataset_path": str(pkl_path),
                "name_dataset": name_dataset,
                "normalize_y": normalize_y,
                "random_state": random_state,
            },
        },
        model_path,
    )
    print(f"💾 Model saved to: {model_path}")

    return {
        "MSE": float(mse),
        "RMSE": rmse,
        "NRMSE": nrmse,
        "MAE": float(mae),
        "R2": float(r2),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "model_path": str(model_path),
        "kernel_": str(gpr.kernel_),
    }