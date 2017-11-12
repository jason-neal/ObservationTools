"""Test code for visibility.py."""

import pytest

from observationtools.visibility import decdeg2dms


@pytest.mark.parametrize("test_input,expected", [
    (0, (0, 0, 0)),
    (33.33, (33, 19, 48)),
    (-45.5125, (-45, 30, 45)),
    (60.60, (60, 36, 0)),
    (-90, (-90, 0, 0)),
])
def test_decdeg2dms(test_input, expected):
    """Test of converting decimal declination."""
    # Possibly limit to +/-90 deg somewhere.
    assert decdeg2dms(test_input) == expected
