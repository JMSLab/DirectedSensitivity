import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import specification as spec

RESULTS_PATH = config.output('estimate.json')


def Main():
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    baseline = Estimate()
    with open(RESULTS_PATH, 'w') as stream:
        json.dump(baseline, stream, indent=2)
        stream.write('\n')
    print(f'wrote {RESULTS_PATH}')
    for name, value in baseline.items():
        print(f'{name:22s} {value!r}')


def Estimate():
    primitive = spec.PUBLISHED_MODEL
    shape = {charity: {'mu': primitive[charity].mu, 'sigma': primitive[charity].sigma}
             for charity in spec.CHARITIES}
    return {**shape,
            'criterion':  spec.Criterion(primitive),
            'target_lar': spec.TargetLaRabida(primitive),
            'target_ecu': spec.TargetECU(primitive)}


if __name__ == '__main__':
    Main()
