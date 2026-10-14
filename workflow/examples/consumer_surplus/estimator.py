"""The criterion, and the least-squares fit of the isoelastic working model."""
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from consumer_surplus import PRICE_RANGE

HERE = Path(__file__).resolve().parent
DATA = HERE / 'demand.csv'
ESTIMATES = HERE / 'estimates.csv'
TOLERANCE = 1e-12
SLOPE_FLOOR = 1e-8


def Main():
    prices, shares = LoadDemand()
    intercept, slope = FitIsoelastic(prices, shares)
    with open(ESTIMATES, 'w', newline='') as stream:
        csv.writer(stream).writerows([('intercept', 'slope'),
                                      (f'{intercept:.10g}', f'{slope:.10g}')])


def LoadDemand(path=DATA):
    data = np.genfromtxt(path, delimiter=',', names=True)
    return data['price'], data['share']


def Criterion(demand, prices, shares):
    """The mean squared gap between the observed shares and demand at the observed prices."""
    return float(np.mean((shares - demand(prices)) ** 2))


def Isoelastic(intercept, slope):
    return lambda price: np.exp(intercept - slope * np.log(price))


def FitIsoelastic(prices, shares):
    """Minimize the criterion over slope > 0 and intercept <= slope log P0."""
    slope, intercept = np.polyfit(np.log(prices), np.log(shares), 1)
    result = minimize(lambda x: Criterion(Isoelastic(*x), prices, shares),
                      x0=[intercept, -slope], method='SLSQP',
                      bounds=[(None, None), (SLOPE_FLOOR, None)],
                      constraints=[{'type': 'ineq',
                                    'fun': lambda x: x[1] * np.log(PRICE_RANGE[0]) - x[0]}],
                      options={'ftol': TOLERANCE, 'maxiter': 1000})
    if not result.success:
        raise RuntimeError(f'the fit did not converge: {result.message}')
    return result.x


if __name__ == '__main__':
    Main()
