"""Checks a proposed basis class should pass before its bounds are believed."""


def NestsBaseline(problem, baseline, tolerance=1e-8):
    """Whether ``baseline`` satisfies the constraints and reproduces the estimate."""
    criterion_gap = _Gap(problem.criterion(baseline), problem.criterion_at_estimate)
    target_gap = _Gap(problem.target(baseline), problem.target_at_estimate)
    violated = problem.ViolatedConstraints(baseline)
    return {'nests': bool(not violated and criterion_gap <= tolerance
                          and target_gap <= tolerance),
            'feasible': problem.IsFeasible(baseline),
            'criterion_gap': criterion_gap,
            'target_gap': target_gap,
            'violated': violated}


def _Gap(value, reference):
    return abs(float(value) - reference) / max(1.0, abs(reference))
