from itertools import product
from pathlib import Path
import pickle

import joblib
import numpy as np



def compute_relative_densities_dataset(lattice_cell, step_radius: float = 0.01, range_radius: tuple = (0.00, 0.1),
                                       name_dataset: str = None):
    """
    Compute relative densities for a range of radii in the lattice.

    Parameters:
    -----------
    lattice : Lattice
        The lattice object containing the cell structure.
    step_radius : float
        Step size for radius values.
    range_radius : tuple
        Tuple specifying the (min, max) radius values.
    """
    if lattice_cell.get_number_cells() != 1:
        raise ValueError("The lattice must contain exactly one cell.")

    dataset_already_exists = False
    path_dataset = Path(__file__).parents[2] / "data" / "outputs" / "relative_densities" / "data"

    if name_dataset is None:
        geom_cell = lattice_cell.cells[0].geom_types
        name_dataset = "RelativeDensities"
        for geom in geom_cell:
            name_dataset += f"_{geom}"

    if (path_dataset / name_dataset).exists():
        print(f"The dataset {name_dataset} already exists in {path_dataset}.")
        dataset_already_exists = True

    if not dataset_already_exists:
        min_radius, max_radius = range_radius
        radius_values = np.arange(min_radius, max_radius + step_radius, step_radius)

        n_geom = len(lattice_cell.cells[0].geom_types)
        all_combinations = list(product(radius_values, repeat=n_geom))

        relative_densities = {}
        for combo in all_combinations:
            if sum(combo) > 0.001:
                print(f"Computing for radii: {combo}")
                lattice_cell.change_beam_radius(list(combo))
                relative_density = lattice_cell.generate_mesh_lattice_Gmsh(volume_computation=True,
                                                                           cut_mesh_at_boundary=True,
                                                                           save_STL=False,
                                                                           only_relative_density=True)
                relative_densities[combo] = relative_density
                print(f"Relative density: {relative_density}")
        save_dataset(path_dataset, name_dataset, relative_densities)
    else:
        print("Dataset already exists.")

def save_dataset(path_dataset, name_dataset, relative_densities_dict):
    """Save the dataset as a pickle file."""
    path_dataset.mkdir(parents=True, exist_ok=True)
    with open(path_dataset / f"{name_dataset}.pkl", "wb") as f:
        pickle.dump(relative_densities_dict, f)

    print(f"Dataset saved as pickle at {path_dataset / (name_dataset + '.pkl')}")

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


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error


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