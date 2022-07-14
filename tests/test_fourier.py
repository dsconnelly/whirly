import numpy as np
import pytest

from whirly.fourier import FourierField
from whirly.utils import make_grid

def make_func(k=1, ell=2):
    def func(x, y):
        return np.sin(k * 2 * np.pi * x) + 0.5 * np.cos(ell * 2 * np.pi * y)

    return func

def test_add_fields():
    a = FourierField.from_func(make_func(), 101)
    b = FourierField.from_func(make_func(2, 1), 101)

    assert np.allclose(a.real + b.real, (a + b).real)

def test_real():
    x, y = make_grid(101, 1)
    real = make_func()(x, y)
    field = FourierField.from_real(real)

    assert np.allclose(real, field.real)

def test_m_must_be_odd():
    hats = np.random.rand(100, 100)
    with pytest.raises(ValueError):
        _ = FourierField(hats, 1)

def test_multiply_fields():
    a = FourierField.from_func(make_func(), 101)
    b = FourierField.from_func(make_func(), 101)

    assert np.allclose(a.real * b.real, (a * b).real)

def test_resample_coarser():
    a = FourierField.from_func(make_func(), 101).resample(75)
    b = FourierField.from_func(make_func(), 75)

    assert np.allclose(a.hats, b.hats)

def test_resample_finer():
    a = FourierField.from_func(make_func(), 101).resample(151)
    b = FourierField.from_func(make_func(), 151)

    assert np.allclose(a.hats, b.hats)
