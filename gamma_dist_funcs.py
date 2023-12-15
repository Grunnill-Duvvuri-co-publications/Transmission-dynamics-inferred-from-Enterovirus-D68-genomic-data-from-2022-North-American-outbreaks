"""Functions for creating dataframes of gamma distribution and making gama distrbution related calcualations."""

from scipy.stats import gamma 
import numpy as np
import pandas as pd


def gamma_from_mean(mu, sd=None, variance=None):
    """Convert mean to shape (a) and scale of gama distribution.

    See:
        Bolker, Benjamin M. 2008. “Gamma.” In Ecological Models and Data in R, 131–133. Princeton University Press.
        https://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
        https://stats.stackexchange.com/questions/342639/how-to-find-alpha-and-beta-from-a-gamma-distribution
        

    Parameters
    ----------
    mean : _type_
        _description_
    sd : _type_, optional
        _description_, by default None
    variance : _type_, optional
        _description_, by default None

    Returns
    -------
    _type_
        _description_
    """
    if variance is None:
        if sd is None:
            variance =1
        else:
            variance = sd**2
    return {'a':(mu**2)/variance, 'scale':variance/mu}

def mean_rate_years_gamma_df(mu, min_range=1, max_range=500, step=1, sd=None, variance=None, return_gamma_dict=False):
    """Generate dataframe outlining pdf of gamma distrubtion for mean yearly rate.

    Parameters
    ----------
    mu : float or int
        Mean in years
    min_range : int, optional
        Min range value in years, by default 1
    max_range : int, optional
        Max range value in years, by default 500
    sd : float or int, optional
        Standard deviation in years., by default None
    variance : float or int, optional
        Variance in years, by default None

    Returns
    -------
    pandas.DataFrame
        Dataframe outlining pdf of gamma distrubtion for mean yearly rate.
    """
    rate_years = np.arange(min_range, max_range+1, step=step)
    location_0s = np.where(rate_years==0)
    rate_days = rate_years/365
    gamma_dict = gamma_from_mean(mu=mu,sd=sd, variance=variance)
    probabilities = gamma.pdf(x=rate_years,**gamma_dict)
    probabilities[location_0s] = 0
    if return_gamma_dict:
        return gamma_dict, pd.DataFrame({'Rate per Year':rate_years, 'Rate per Day':rate_days, 'Probability': probabilities})
    else:
        pd.DataFrame({'Rate per Year':rate_years, 'Rate per Day':rate_days, 'Probability': probabilities})