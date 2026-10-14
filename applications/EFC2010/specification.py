import json
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

MODEL = config.CONFIG['model']
GRIDS = config.CONFIG['grids']
TOL   = config.CONFIG['tolerances']

PRICE_INTERVAL = tuple(MODEL['price_interval'])
TARGET_GRID    = np.linspace(*PRICE_INTERVAL, GRIDS['target_points'])


class LinearDemand:
    def __init__(self, intercept, slope):
        self.intercept, self.slope = float(intercept), float(slope)

    def Value(self, price):
        return self.intercept + self.slope * np.asarray(price, dtype=float)

    def Derivative(self, price):
        return np.full(np.shape(np.asarray(price, dtype=float)), self.slope)


def LoadData(path=config.PATHS['data']):
    table = pd.read_csv(path, thousands=',')
    employees = table['number_of_employees'].values.astype(float)
    share = table['fraction_chose_contract_H'].values.astype(float)
    return {'p':  table['relative_price'].values.astype(float),
            'D':  share,
            'c':  table['average_incremental_cost_contract_H'].values.astype(float),
            'w1': employees,
            'w2': employees * share}


def FitLine(x, y, weights):
    return tuple(sm.WLS(y, sm.add_constant(x), weights=weights).fit().params)


DATA = LoadData()
PRICE, SHARE, COST, EMPLOYEES = DATA['p'], DATA['D'], DATA['c'], DATA['w1']
WEIGHTS = MODEL['criterion_scale'] * EMPLOYEES / EMPLOYEES.sum()

AVG_COST_INTERCEPT, AVG_COST_SLOPE = FitLine(PRICE, COST, DATA['w2'])
COMPETITIVE_PRICE = AVG_COST_INTERCEPT / (1.0 - AVG_COST_SLOPE)


def Main():
    baseline = Baseline()
    print(f'demand intercept, slope    {baseline.intercept!r}, {baseline.slope!r}')
    print(f'avg cost intercept, slope  {AVG_COST_INTERCEPT!r}, {AVG_COST_SLOPE!r}')
    print(f'competitive price          {COMPETITIVE_PRICE!r}')
    print(f'criterion                  {Criterion(baseline)!r}')
    print(f'efficient price            {EfficientPrice(baseline)!r}')
    print(f'target                     {Target(baseline)!r}')


def Baseline(path=None):
    path = config.output('estimate.json') if path is None else path
    with open(path) as stream:
        estimate = json.load(stream)
    return LinearDemand(estimate['demand_intercept'], estimate['demand_slope'])


def Criterion(primitive):
    residual = SHARE - primitive.Value(PRICE)
    return float(WEIGHTS @ residual ** 2)


def Target(primitive):
    slope = np.asarray(primitive.Derivative(TARGET_GRID), dtype=float)
    if np.any(slope == 0.0) or not np.all(np.isfinite(slope)):
        raise ValueError('the derivative vanishes or is not finite on '
                         f'{PRICE_INTERVAL}, so marginal cost is undefined')

    value, _ = quad(lambda price: WelfareIntegrand(primitive, price),
                    COMPETITIVE_PRICE, EfficientPrice(primitive),
                    epsabs=TOL['integration'])
    return value


def WelfareIntegrand(primitive, price):
    return float((1.0 - AVG_COST_SLOPE) * (price - COMPETITIVE_PRICE)
                 * primitive.Derivative(price)
                 - AVG_COST_SLOPE * primitive.Value(price))


def AverageCost(price):
    return AVG_COST_INTERCEPT + AVG_COST_SLOPE * np.asarray(price, dtype=float)


def MarginalCost(primitive, price):
    return (AverageCost(price)
            + AVG_COST_SLOPE * primitive.Value(price) / primitive.Derivative(price))


def EfficientPrice(primitive):
    def gap(price):
        return float(price - MarginalCost(primitive, price))

    if gap(0.0) >= 0.0:
        return 0.0
    if gap(COMPETITIVE_PRICE) <= 0.0:
        return COMPETITIVE_PRICE
    return brentq(gap, 0.0, COMPETITIVE_PRICE, xtol=TOL['root'])


if __name__ == '__main__':
    Main()
