"""The object of interest: the loss in consumer surplus from a price increase."""
from scipy.integrate import quad

PRICE_RANGE = (100.0, 200.0)
TOLERANCE = 1e-10


def ConsumerSurplusLoss(demand, low=PRICE_RANGE[0], high=PRICE_RANGE[1]):
    """The area under the demand curve between the two prices."""
    value, _ = quad(demand, low, high, epsabs=TOLERANCE, epsrel=TOLERANCE)
    return float(value)
