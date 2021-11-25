from pdiffutils import DiffractionDataPoint
import math


# remember, I've overridden ==
def is_all_equal(p1, p2):
    return bool(
        math.isclose(p1.x, p2.x)
        and math.isclose(p1.y, p2.y)
        and math.isclose(p1.e, p2.e)
    )


def test_construction():
    ddp1 = DiffractionDataPoint(5.01, 25.0)
    ddp2 = DiffractionDataPoint(5.01, -25.0)
    ddp3 = DiffractionDataPoint(5.01, 25.0, 3.6)
    ddp4 = DiffractionDataPoint(5.01, 25.0, -3.6)

    assert ddp1.x == 5.01
    assert ddp1.y == 25.0
    assert ddp1.e == 5.0

    assert ddp2.x == 5.01
    assert ddp2.y == -25.0
    assert ddp2.e == 5.0

    assert ddp3.x == 5.01
    assert ddp3.y == 25.0
    assert ddp3.e == 3.6

    assert ddp4.x == 5.01
    assert ddp4.y == 25.0
    assert ddp4.e == 3.6


def test_negate():
    ddp1 = DiffractionDataPoint(5.01, 25.0)
    ddp1.negate()
    assert is_all_equal(ddp1, DiffractionDataPoint(-5.01, 25.0))


def test_zeroOffset():
    ddp1 = DiffractionDataPoint(5.01, 25.0)
    ddp1.zeroOffset(-0.01)
    assert is_all_equal(ddp1, DiffractionDataPoint(5.00, 25.0))


def test_average_with():
    ddp1 = DiffractionDataPoint(5.01, 3.0)
    ddp2 = DiffractionDataPoint(5.01, 4.0)
    ddp3 = DiffractionDataPoint(5.01, 5.0)
    ddp4 = DiffractionDataPoint(5.01, 7.0)

    ddp1.average_with([ddp2, ddp3, ddp4])

    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, 4.75, math.sqrt(19.0) / 4))


def test_comparisons():
    ddp1 = DiffractionDataPoint(5.01, 25.0)
    ddp2 = DiffractionDataPoint(5.01, 27.0)
    ddp3 = DiffractionDataPoint(5.02, 25.0)
    ddp4 = DiffractionDataPoint(5.02, 29.0)

    assert ddp1 < ddp3
    assert ddp1 <= ddp2
    assert ddp4 > ddp1
    assert ddp3 >= ddp2
    assert ddp1 == ddp2
    assert ddp1 != ddp3


def test_operators():
    ddp1 = DiffractionDataPoint(5.01, 9.0, 3.0)
    ddp2 = DiffractionDataPoint(5.01, 16.0, 4.0)

    assert is_all_equal(ddp1 * 2, DiffractionDataPoint(5.01, 18.0, 6.0))
    assert is_all_equal(ddp1 * -2, DiffractionDataPoint(5.01, -18.0, 6.0))
    assert is_all_equal(ddp1 / 2, DiffractionDataPoint(5.01, 4.5, 1.5))
    assert is_all_equal(ddp1 / -2, DiffractionDataPoint(5.01, -4.5, 1.5))
    assert is_all_equal(ddp1 // 2, DiffractionDataPoint(5.01, 4, 1))
    assert is_all_equal(ddp1 // -2, DiffractionDataPoint(5.01, -5, 1))
    assert is_all_equal(ddp1 + 2, DiffractionDataPoint(5.01, 11.0, 3.0))
    assert is_all_equal(ddp1 - 2, DiffractionDataPoint(5.01, 7.0, 3.0))
    assert is_all_equal(ddp1 + ddp2, DiffractionDataPoint(5.01, 25.0, 5.0))
    assert is_all_equal(ddp1 - ddp2, DiffractionDataPoint(5.01, -7.0, 5.0))


def test_ioperators():
    ddp1 = DiffractionDataPoint(5.01, 9.0, 3.0)
    ddp2 = DiffractionDataPoint(5.01, 16.0, 4.0)

    ddp1 *= 2
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, 18.0, 6.0))
    ddp1 *= -3
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, -54.0, 18.0))

    ddp1 /= 2
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, -27, 9.0))
    ddp1 /= -3
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, 9.0, 3.0))

    ddp1 //= 2
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, 4, 1))
    ddp2 //= -3
    assert is_all_equal(ddp2, DiffractionDataPoint(5.01, -6, 1))

    ddp1 += 2
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, 6, 1))
    ddp2 -= 2
    assert is_all_equal(ddp2, DiffractionDataPoint(5.01, -8, 1))

    ddp1 += ddp2
    assert is_all_equal(ddp1, DiffractionDataPoint(5.01, -2, math.sqrt(2)))
    ddp2 -= ddp1
    assert is_all_equal(ddp2, DiffractionDataPoint(5.01, -6, math.sqrt(3)))
