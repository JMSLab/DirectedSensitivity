import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import specification as spec

RESULTS_PATH = config.output('replicate.txt')

# Table 2, column (5) and Table 4, column (1).
PUBLISHED = {'present_bias': 0.755, 'cost': 2.480, 'log_likelihood': -66290.55}


def Main():
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    outfile = open(RESULTS_PATH, 'w')
    sys.stdout = _Tee(sys.__stdout__, outfile)

    recovered = Recover()
    PrintData()
    PrintEstimate(recovered)
    tolerance = 0.005
    CheckAssertions(recovered, tolerance)


class _Tee:
    def __init__(self, *files):
        self._files = files

    def write(self, text):
        for f in self._files:
            f.write(text)

    def flush(self):
        for f in self._files:
            f.flush()


def Recover():
    return {'present_bias': float(spec.ESTIMATE[0]), 'cost': -float(spec.ESTIMATE[1]),
            'log_likelihood': -spec.Criterion(spec.WORKING_MODEL)}


def PrintData():
    print(f'Decision days                  {spec.SEASON_DAYS}')
    print(f'Filers at risk on day 0        {spec.AT_RISK[0]:.0f}')
    print(f'Filings over the decision days {spec.FILED.sum():.0f}')
    print(f'Person-days at risk            {spec.AT_RISK.sum():.0f}')


def PrintEstimate(recovered):
    print(f"\n{'':<16}{'recovered':>16}{'published':>16}")
    for key, digits in (('present_bias', 6), ('cost', 6), ('log_likelihood', 4)):
        print(f'{key:<16}{recovered[key]:>16.{digits}f}{PUBLISHED[key]:>16.{digits}f}')
    print(f"\nRounded to the printed precision: "
          f"beta {recovered['present_bias']:.3f}, cost {recovered['cost']:.3f}, "
          f"log-likelihood {recovered['log_likelihood']:.2f}")


def CheckAssertions(recovered, tolerance):
    for key in ('present_bias', 'cost', 'log_likelihood'):
        assert np.isclose(recovered[key], PUBLISHED[key], rtol=tolerance), (
            f'FAIL {key}: computed {recovered[key]:.6g}, '
            f'target {PUBLISHED[key]:.6g} (rtol={tolerance})')


if __name__ == '__main__':
    Main()
