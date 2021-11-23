import numpy as np
import pandas as pd
from scipy import constants, stats

EPSILON = 1e-10


def great_circle_dist1(r1, d1, r2, d2):
    return np.arccos(
        np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    )


def precise_dist(ra1, dec1, ra2, dec2):
    d1 = np.deg2rad(dec1)
    d2 = np.deg2rad(dec2)
    r1 = np.deg2rad(ra1)
    r2 = np.deg2rad(ra2)
    dist = great_circle_dist1(r1, d1, r2, d2)
    dist = dist * 180.0 / np.pi
    return dist


# From `https://gist.github.com/bshishov/5dc237f59f019b26145648e2124ca1c9`
def _error(actual: np.ndarray, predicted: np.ndarray):
    """Simple error"""
    return actual - predicted


def mse(actual: np.ndarray, predicted: np.ndarray):
    """Mean Squared Error"""
    return np.mean(np.square(_error(actual, predicted)))


# def rmse(xx, yy):
#     return mean_squared_error(xx, yy, squared=False) / (np.max(xx) - np.min(xx)) * 100
def rmse(actual: np.ndarray, predicted: np.ndarray):
    """Root Mean Squared Error"""
    return np.sqrt(mse(actual, predicted))


def nrmse(actual, predicted):
    """Normalized Root Mean Squared Error"""
    return rmse(actual, predicted) / (actual.max() - actual.min())


def mape_maape(actual, predicted, avg=True):
    try:
        difference = np.abs(actual - predicted)
        error = (difference / actual) * 100

    except Exception:
        print("using maape")
        error = np.arctan(np.abs((actual - predicted) / (actual + EPSILON))) * 100
    avg_error = np.mean(error)
    if avg:
        return avg_error
    else:
        return avg_error, error


def get_metric(iono_mag, pca):
    return 25 * iono_mag + 64 * pca * (pca - 0.6) if pca > 0.6 else 25 * iono_mag


def scale_offsets(input_offsets, input_frequency=154.235, output_frequency=200):
    input_offsets *= (input_frequency ** 2) / (output_frequency ** 2)
    return input_offsets


def get_iono_type(m, p):
    if m < 0.14 and p < 0.63:
        itype = 1
    elif m > 0.14 and p < 0.7:
        itype = 2
    elif m < 0.14 and p > 0.63:
        itype = 3
    elif m > 0.14 and p > 0.7:
        itype = 4
    else:
        print("weird")
    return itype


def iqr_bounds(arr):
    # Calculating min, max, Q1, and Q3 interquartile ranges
    iqr = stats.iqr(arr)
    q1 = np.quantile(a=arr, q=0.25)
    minimum = q1 - 1.5 * iqr
    q3 = np.quantile(a=arr, q=0.75)
    maximum = q3 + 1.5 * iqr

    if minimum < arr.min():
        min_fence = arr.min()
    else:
        min_fence = minimum
    if maximum > arr.max():
        max_fence = arr.max()
    else:
        max_fence = maximum

    return min_fence, max_fence


def thermal_noise(
    effective_collecting_area=21.4,  # m^2
    dnu=40000.0000000298 * 2,  # RTS uses 80KHz nu channels
    integration_time=8,  # s
    Tsys=316,  # K #300 #266+50
    n_ants=123,
):
    n_baselines = 0.5 * n_ants * (n_ants - 1)
    sigma = (
        2
        * 1e26  # Conversion to Jys
        * constants.k
        * Tsys
        / effective_collecting_area
        / np.sqrt(
            n_baselines * dnu * integration_time
        )  # goes down as 1/sqrt(n_baselines)
    )
    return sigma


def get_2sigma_sources(obsid, df):
    tisis = get_sky_temp(obsid) + 50
    print("Tsys:", tisis)
    thermal_noise_level = thermal_noise(Tsys=tisis)
    two_sigma_df = df[df["flux"] >= thermal_noise_level * 2]

    return two_sigma_df.name.to_list()


def get_sky_temp(obsid):
    # csvpath = "/home/kariuki/gara_pipe/sampled_5obsids_per_iono_type.csv"
    csvpath = "/astro/mwaeor/kchege/pipeline2/dataset_selection/csvs/final_lst0_clean_lowband_info_db.csv"
    df = pd.read_csv(csvpath, usecols=["obsid", "sky_temp"])
    temp = df[df["obsid"] == int(obsid)]["sky_temp"].values[0]

    return float(temp)
