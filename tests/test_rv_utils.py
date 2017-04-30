import pytest
import numpy as np
from utils.rv_utils import mean_anomaly
from utils.rv_utils import true_anomaly

from hypothesis import strategies as st
from hypothesis import given, example, assume


# issue with limits  0-pi only
@given(st.lists(st.floats(min_value=0, max_value=np.pi), min_size=1), st.floats(min_value=0.05, max_value=0.99))
def test_trueanomaly(ma, ecc):

    ma = np.asarray(ma)
    assume(np.all(np.abs(ma) > 0.0001))
    ta = true_anomaly(ma, ecc)
    # cos(E) = (e + cos(ta)) / (1 + e*cos(ta))
    E = np.arccos((ecc + np.cos(ta)) /(1 + ecc * np.cos(ta)))
    # ma = E - e*sin(E)
    print("VALUES", ma, (E - ecc * np.sin(E)))
    assert np.allclose(ma, E - ecc * np.sin(E), rtol=0.05)
    assert len(ta) == len(ma)


# issue with limits 0-pi only
@given(st.floats(min_value=0, max_value=np.pi), st.floats(min_value=0.01, max_value=0.99))
@example(2, 0.5)   # example with an integer
def test_trueanomaly_with_scalar(ma, ecc):
    assume(abs(ma) > 0.001)
    ta = true_anomaly(ma, ecc)
    # Contrast Result to E from ta
    # cos(E) = (e + cos(ta)) / (1 + e*cos(ta))
    E = np.arccos((ecc + np.cos(ta)) / (1 + ecc * np.cos(ta)))
    # ma = E - e*sin(E)
    MA = E - ecc * np.sin(E)
    assert np.allclose(ma, E - ecc * np.sin(E))
    assert len(ta) == 1


@given(st.floats(min_value=0.01, max_value=0.99))
def test_trueanomaly_errors(ecc):

    with pytest.raises(TypeError):
        true_anomaly([], ecc)

    with pytest.raises(ValueError):
        true_anomaly(np.array([]), ecc)


@given(st.lists(st.floats(), min_size=1), st.floats(), st.floats(min_value=0.01))
def test_mean_anomaly(t, t0, p):
    """Mean anomaly is an angle, doen't have a constraint value."""
    t = np.array(t)
    ma = mean_anomaly(t, t0, p)

    assert len(t) == len(ma)
    assert isinstance(t, np.ndarray)