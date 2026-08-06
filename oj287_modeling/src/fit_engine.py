import numpy as np


def chi_squared(observed, predicted, sigma):
    observed = np.asarray(observed)
    predicted = np.asarray(predicted)
    sigma = np.asarray(sigma)
    return np.sum(((observed - predicted) / sigma) ** 2)


def log_likelihood(observed, predicted, sigma):
    observed = np.asarray(observed)
    predicted = np.asarray(predicted)
    sigma = np.asarray(sigma)
    chi2 = ((observed - predicted) / sigma) ** 2
    return -0.5 * (np.sum(chi2) + np.sum(np.log(2 * np.pi * sigma ** 2)))
