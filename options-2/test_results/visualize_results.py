import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from os.path import join


# without bias
def apply_OLS(X: np.ndarray, y: np.ndarray):
    X_t = X.T
    theta = np.dot(np.linalg.inv(np.dot(X_t, X)), np.dot(X_t, y))
    mse = np.mean(np.square(y - np.dot(X, theta)))
    return theta, mse


def plot_results():
    results = pd.read_csv("n_vs_norm.csv")

    X = np.log(results['n'].values.reshape(-1, 1))
    X = np.c_[np.ones(len(X)), X]

    y_ge = np.log(results['norm_ge'].values.reshape(-1, 1))
    theta_ge, _ = apply_OLS(X, y_ge)

    y_gsl = np.log(results['norm_gsl'].values.reshape(-1, 1))
    theta_gsl, _ = apply_OLS(X, y_gsl)

    fig, ax = plt.subplots()

    ax.plot(X[:, 1], y_ge, label='GE')
    ax.plot(X[:, 1], np.dot(X, theta_ge), label=f'slope = {theta_ge[1][0]:.5f} (GE)')
    # ax.plot(X[:, 1], y_gsl, label='GSL')
    # ax.plot(X[:, 1], np.dot(X, theta_gsl), label=f'slope = {theta_gsl[1][0]:.5f} (GSL)')

    ax.set(title="error(N) plot", xlabel='ln(N)', ylabel='ln(error(N))')
    ax.legend()

    fig.tight_layout()
    fig.savefig(fname="n_vs_norm.png")


plot_results()
