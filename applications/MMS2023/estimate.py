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
    present_bias, location = (float(value) for value in spec.ESTIMATE)
    return {'present_bias': present_bias,
            'location':     location,
            'criterion':    spec.Criterion(spec.WORKING_MODEL),
            'target':       spec.Target(spec.WORKING_MODEL)}


if __name__ == '__main__':
    Main()
