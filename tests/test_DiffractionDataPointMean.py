from pdiffutils import DiffractionDataPoint, DiffractionDataPointMean
import math
import pytest


# remember, I've overridden ==
def is_all_equal(p1, p2):
    return bool(
        math.isclose(p1.x, p2.x)
        and math.isclose(p1.y, p2.y)
        and math.isclose(p1.e, p2.e)
    )


def test_construction():
    ddpmean = DiffractionDataPointMean()
    assert ddpmean.n == 0
    assert math.isnan(ddpmean.x)
    assert math.isnan(ddpmean.m1_y)
    assert math.isnan(ddpmean.dev_y)
    assert math.isnan(ddpmean.nDev_y)
    assert math.isnan(ddpmean.m1_e)
    assert math.isnan(ddpmean.dev_e)
    assert math.isnan(ddpmean.nDev_e)


def test_increment():
    ddp1 = DiffractionDataPoint(5.01, 3.0)
    ddp2 = DiffractionDataPoint(5.01, 4.0)
    ddp3 = DiffractionDataPoint(5.01, 5.0)
    ddp4 = DiffractionDataPoint(5.02, 7.0)
    ddpmean = DiffractionDataPointMean()

    ddpmean.increment(ddp1)
    ddpmean.increment(ddp2)
    ddpmean.increment(ddp3)

    assert ddpmean.n == 3
    assert math.isclose(ddpmean.m1_y, 4.0)
    assert math.isclose(ddpmean.m1_e, 4.0)

    with pytest.raises(ValueError) as e_info:
        ddpmean.increment(ddp4)

def test_get_results():
    ddp1 = DiffractionDataPoint(5.01, 3.0)
    ddp2 = DiffractionDataPoint(5.01, 4.0)
    ddp3 = DiffractionDataPoint(5.01, 5.0)
    ddpmean = DiffractionDataPointMean()

    ddpmean.increment(ddp1)
    ddpmean.increment(ddp2)
    ddpmean.increment(ddp3)

    assert is_all_equal(ddpmean.get_result(), DiffractionDataPoint(5.01, (3.0+4.0+5.0)/3, math.sqrt(3.0+4.0+5.0)/3))


def test_clear():
    ddp1 = DiffractionDataPoint(5.01, 3.0)
    ddp2 = DiffractionDataPoint(5.01, 4.0)
    ddp3 = DiffractionDataPoint(5.01, 5.0)
    ddpmean = DiffractionDataPointMean()
    ddpmean.increment(ddp1)
    ddpmean.increment(ddp2)
    ddpmean.increment(ddp3)

    ddpmean.clear()
    assert ddpmean.n == 0
    assert math.isnan(ddpmean.x)
    assert math.isnan(ddpmean.m1_y)
    assert math.isnan(ddpmean.dev_y)
    assert math.isnan(ddpmean.nDev_y)
    assert math.isnan(ddpmean.m1_e)
    assert math.isnan(ddpmean.dev_e)
    assert math.isnan(ddpmean.nDev_e)

