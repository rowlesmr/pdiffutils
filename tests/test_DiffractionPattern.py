import copy

import pytest

from pdiffutils import DiffractionPattern, DiffractionDataPoint
import math
import operator


# remember, I've overridden ==
def is_all_equal(d1, d2):
    return not any(
        not math.isclose(p1.x, p2.x)
        or not math.isclose(p1.y, p2.y)
        or not math.isclose(p1.e, p2.e)
        for p1, p2 in zip(d1.diffpat, d2.diffpat)
    )


def make_dp(angle: float = 0, intensity: float = 0, error: float = 0, op: operator = None) -> DiffractionPattern:
    if op is None:
        angle = 0
        intensity = 0
        error = 0
        op = operator.add

    diffpat2 = [DiffractionDataPoint(op(5.00, angle), op(4.1, intensity), op(math.sqrt(4.1), error)),
                DiffractionDataPoint(op(5.01, angle), op(2.1, intensity), op(math.sqrt(2.1), error)),
                DiffractionDataPoint(op(5.02, angle), op(3.1, intensity), op(math.sqrt(3.1), error)),
                DiffractionDataPoint(op(5.03, angle), op(4.1, intensity), op(math.sqrt(4.1), error)),
                DiffractionDataPoint(op(5.04, angle), op(6.1, intensity), op(math.sqrt(6.1), error)),
                DiffractionDataPoint(op(5.05, angle), op(5.1, intensity), op(math.sqrt(5.1), error))]
    return DiffractionPattern(diffpat=diffpat2)


def test_construction():
    diffpat = [DiffractionDataPoint(5.00, 2.1),
               DiffractionDataPoint(5.01, 4.1),
               DiffractionDataPoint(5.02, 4.1),
               DiffractionDataPoint(5.03, 3.1),
               DiffractionDataPoint(5.04, 6.1)]

    dp = DiffractionPattern(diffpat=diffpat)

    assert dp.xs == [5.0, 5.01, 5.02, 5.03, 5.04]
    assert dp.ys == [2.1, 4.1, 4.1, 3.1, 6.1]
    assert dp.es == [math.sqrt(2.1), math.sqrt(4.1), math.sqrt(4.1), math.sqrt(3.1), math.sqrt(6.1)]
    assert math.isclose(dp.ave_step_size, 0.01)


def test_negate():
    dp1 = make_dp()
    dp2 = make_dp(-1, 1, 1, operator.mul)
    dp1.negate()
    assert is_all_equal(dp1, dp2)


def test_reverse():
    diffpat1 = [DiffractionDataPoint(5.00, 4.1),
                DiffractionDataPoint(5.01, 2.1),
                DiffractionDataPoint(5.02, 3.1),
                DiffractionDataPoint(5.03, 4.1),
                DiffractionDataPoint(5.04, 6.1),
                DiffractionDataPoint(5.05, 5.1)]
    diffpat2 = [DiffractionDataPoint(5.05, 5.1),
                DiffractionDataPoint(5.04, 6.1),
                DiffractionDataPoint(5.03, 4.1),
                DiffractionDataPoint(5.02, 3.1),
                DiffractionDataPoint(5.01, 2.1),
                DiffractionDataPoint(5.00, 4.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    dp2 = DiffractionPattern(diffpat=diffpat2)
    dp1.reverse()
    assert is_all_equal(dp1, dp2)


def test_zeroOffset():
    dp1 = make_dp()
    dp2 = make_dp(op=operator.add, angle=0.1)
    dp1.zero_offset(0.1)
    assert is_all_equal(dp1, dp2)

    diffpat2 = [DiffractionDataPoint(5.00, 4.1),
                DiffractionDataPoint(5.11, 2.1),
                DiffractionDataPoint(5.22, 3.1),
                DiffractionDataPoint(5.33, 4.1),
                DiffractionDataPoint(5.44, 6.1),
                DiffractionDataPoint(5.55, 5.1)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    dp1.zero_offset([-0.1, 0, .1, .2, .3, .4])
    assert is_all_equal(dp1, dp2)

    with pytest.raises(ValueError) as e_info:
        dp1.zero_offset([-0.1, 0, .1, .2])


def test_len():
    dp1 = make_dp()
    assert len(dp1) == 6


def test_trim():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1)]
    dp2 = DiffractionPattern(diffpat=diffpat2)

    assert is_all_equal(dp1.trim(in_place=False), dp1)
    assert is_all_equal(dp1.trim(min_x=5.01, max_x=5.03, in_place=False), dp2)


def test_sort():
    diffpat1 = [DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp2 = DiffractionPattern(diffpat=diffpat2)

    assert is_all_equal(dp1.sort(in_place=False), dp2)

    dp1.sort(in_place=True)
    assert is_all_equal(dp1, dp2)


def test_downsample():
    diffpat1 = [DiffractionDataPoint(5.00, 4.1),
                DiffractionDataPoint(5.01, 2.1),
                DiffractionDataPoint(5.02, 3.1),
                DiffractionDataPoint(5.03, 4.1),
                DiffractionDataPoint(5.04, 6.1),
                DiffractionDataPoint(5.05, 5.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.01, (4.1 + 2.1 + 3.1) / 3, math.sqrt(4.1 + 2.1 + 3.1) / 3),
                DiffractionDataPoint(5.04, (4.1 + 6.1 + 5.1) / 3, math.sqrt(4.1 + 6.1 + 5.1) / 3)]
    dp2 = DiffractionPattern(diffpat=diffpat2)

    assert is_all_equal(dp1.downsample(3, in_place=False), dp2)

    dp1.downsample(3, in_place=True)
    assert is_all_equal(dp1, dp2)


def test_split_on_zero():
    diffpat1 = [DiffractionDataPoint(5.00, 4.1),
                DiffractionDataPoint(5.01, 2.1),
                DiffractionDataPoint(5.02, 3.1),
                DiffractionDataPoint(5.03, 4.1),
                DiffractionDataPoint(5.04, 6.1),
                DiffractionDataPoint(5.05, 5.1)]

    diffpat2 = [DiffractionDataPoint(-5.05, 5.1),
                DiffractionDataPoint(-5.04, 6.1),
                DiffractionDataPoint(-5.03, 4.1),
                DiffractionDataPoint(-5.02, 3.1),
                DiffractionDataPoint(-5.01, 2.1),
                DiffractionDataPoint(-5.00, 4.1)]

    diffpat = copy.deepcopy(diffpat2 + diffpat1)
    dp0 = DiffractionPattern(diffpat=diffpat)
    dp1 = DiffractionPattern(diffpat=diffpat1)
    dp2 = DiffractionPattern(diffpat=diffpat2)
    dp2.negate()
    dp2.reverse()

    dpp, dpn = dp0.split_on_zero()
    assert is_all_equal(dpp, dp1)
    assert is_all_equal(dpn, dp2)  # the negative pattern is returned as the equivalent positive pattern


def test_average_with():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    diffpat2 = [DiffractionDataPoint(5.00, 3.2),
                DiffractionDataPoint(5.01, 4.2),
                DiffractionDataPoint(5.02, 5.2),
                DiffractionDataPoint(5.03, 5.2),
                DiffractionDataPoint(5.04, 7.2)]
    diffpat3 = [DiffractionDataPoint(5.00, 5.3),
                DiffractionDataPoint(5.01, 6.3),
                DiffractionDataPoint(5.02, 5.3),
                DiffractionDataPoint(5.03, 3.3),
                DiffractionDataPoint(5.04, 7.3)]
    diffpat4 = [DiffractionDataPoint(5.00, (2.1 + 3.2 + 5.3) / 3, math.sqrt(2.1 + 3.2 + 5.3) / 3),
                DiffractionDataPoint(5.01, (4.1 + 4.2 + 6.3) / 3, math.sqrt(4.1 + 4.2 + 6.3) / 3),
                DiffractionDataPoint(5.02, (4.1 + 5.2 + 5.3) / 3, math.sqrt(4.1 + 5.2 + 5.3) / 3),
                DiffractionDataPoint(5.03, (3.1 + 5.2 + 3.3) / 3, math.sqrt(3.1 + 5.2 + 3.3) / 3),
                DiffractionDataPoint(5.04, (6.1 + 7.2 + 7.3) / 3, math.sqrt(6.1 + 7.2 + 7.3) / 3)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    dp2 = DiffractionPattern(diffpat=diffpat2)
    dp3 = DiffractionPattern(diffpat=diffpat3)
    dp4 = DiffractionPattern(diffpat=diffpat4)

    assert is_all_equal(dp1.average_with([dp2, dp3], in_place=False), dp4)

    dp1.average_with([dp2, dp3], in_place=True)
    assert is_all_equal(dp1, dp4)


def test_operators():
    dp1 = make_dp()

    assert is_all_equal(dp1 * 2, make_dp(1, 2, 2, operator.mul))
    assert is_all_equal(dp1 * -2, make_dp(1, -2, -2, operator.mul))
    assert is_all_equal(dp1 / 2, make_dp(1, 2, 2, operator.truediv))
    assert is_all_equal(dp1 / -2, make_dp(1, -2, -2, operator.truediv))
    # assert is_all_equal(dp1 // 2, make_dp(1, 2, 2, operator.floordiv))
    # assert is_all_equal(dp1 // -2, make_dp(1, -2, abs(-2), operator.floordiv))
    assert is_all_equal(dp1 + 2, make_dp(0, 2, 0, operator.add))
    assert is_all_equal(dp1 - 2, make_dp(0, 2, 0, operator.sub))

    diffpat3 = [DiffractionDataPoint(5.00, 4.1 + 1, math.sqrt(4.1)),
                DiffractionDataPoint(5.01, 2.1 + 2, math.sqrt(2.1)),
                DiffractionDataPoint(5.02, 3.1 + 3, math.sqrt(3.1)),
                DiffractionDataPoint(5.03, 4.1 + 4, math.sqrt(4.1)),
                DiffractionDataPoint(5.04, 6.1 + 5, math.sqrt(6.1)),
                DiffractionDataPoint(5.05, 5.1 + 6, math.sqrt(5.1))]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    assert is_all_equal(dp1 + [1, 2, 3, 4, 5, 6], dp3)
    assert is_all_equal(dp1 - [-1, -2, -3, -4, -5, -6], dp3)

    dp2 = make_dp()
    diffpat3 = [DiffractionDataPoint(5.00, 4.1 + 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.01, 2.1 + 2.1, math.sqrt(2.1 + 2.1)),
                DiffractionDataPoint(5.02, 3.1 + 3.1, math.sqrt(3.1 + 3.1)),
                DiffractionDataPoint(5.03, 4.1 + 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.04, 6.1 + 6.1, math.sqrt(6.1 + 6.1)),
                DiffractionDataPoint(5.05, 5.1 + 5.1, math.sqrt(5.1 + 5.1))]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    assert is_all_equal(dp1 + dp2, dp3)

    diffpat3 = [DiffractionDataPoint(5.00, 4.1 - 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.01, 2.1 - 2.1, math.sqrt(2.1 + 2.1)),
                DiffractionDataPoint(5.02, 3.1 - 3.1, math.sqrt(3.1 + 3.1)),
                DiffractionDataPoint(5.03, 4.1 - 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.04, 6.1 - 6.1, math.sqrt(6.1 + 6.1)),
                DiffractionDataPoint(5.05, 5.1 - 5.1, math.sqrt(5.1 + 5.1))]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    assert is_all_equal(dp1 - dp2, dp3)

    with pytest.raises(TypeError) as e_info:
        _ = dp1 * dp2
        _ = dp1 / dp2
        _ = dp1 // dp2
        _ = dp1 + (1, 2, 3)
        _ = dp1 - (1, 2, 3)

    with pytest.raises(ValueError) as e_info:
        _ = dp1 + [1, 2, 3]
        _ = dp1 - [1, 2, 3]


def test_ioperators():
    dp1 = make_dp()
    dp1 *= 2
    assert is_all_equal(dp1, make_dp(1, 2, 2, operator.mul))
    dp1 = make_dp()
    dp1 *= -2
    assert is_all_equal(dp1, make_dp(1, -2, -2, operator.mul))
    dp1 = make_dp()
    dp1 /= 2
    assert is_all_equal(dp1, make_dp(1, 2, 2, operator.truediv))
    dp1 = make_dp()
    dp1 /= -2
    assert is_all_equal(dp1, make_dp(1, -2, -2, operator.truediv))
    dp1 = make_dp()
    dp1 //= 2
    # assert is_all_equal(dp1, make_dp(1, 2, 2, operator.floordiv))
    dp1 = make_dp()
    dp1 //= -2
    # assert is_all_equal(dp1, make_dp(1, -2, abs(-2), operator.floordiv))
    dp1 = make_dp()
    dp1 += 2
    assert is_all_equal(dp1, make_dp(0, 2, 0, operator.add))
    dp1 = make_dp()
    dp1 -= 2
    assert is_all_equal(dp1, make_dp(0, 2, 0, operator.sub))

    diffpat3 = [DiffractionDataPoint(5.00, 4.1 + 1, math.sqrt(4.1)),
                DiffractionDataPoint(5.01, 2.1 + 2, math.sqrt(2.1)),
                DiffractionDataPoint(5.02, 3.1 + 3, math.sqrt(3.1)),
                DiffractionDataPoint(5.03, 4.1 + 4, math.sqrt(4.1)),
                DiffractionDataPoint(5.04, 6.1 + 5, math.sqrt(6.1)),
                DiffractionDataPoint(5.05, 5.1 + 6, math.sqrt(5.1))]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    dp1 = make_dp()
    dp1 += [1, 2, 3, 4, 5, 6]
    assert is_all_equal(dp1, dp3)
    dp1 = make_dp()
    dp1 -= [-1, -2, -3, -4, -5, -6]
    assert is_all_equal(dp1, dp3)

    diffpat3 = [DiffractionDataPoint(5.00, 4.1 + 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.01, 2.1 + 2.1, math.sqrt(2.1 + 2.1)),
                DiffractionDataPoint(5.02, 3.1 + 3.1, math.sqrt(3.1 + 3.1)),
                DiffractionDataPoint(5.03, 4.1 + 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.04, 6.1 + 6.1, math.sqrt(6.1 + 6.1)),
                DiffractionDataPoint(5.05, 5.1 + 5.1, math.sqrt(5.1 + 5.1))]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    dp2 = make_dp()

    dp1 = make_dp()
    dp1 += dp2
    assert is_all_equal(dp1, dp3)
    diffpat3 = [DiffractionDataPoint(5.00, 4.1 - 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.01, 2.1 - 2.1, math.sqrt(2.1 + 2.1)),
                DiffractionDataPoint(5.02, 3.1 - 3.1, math.sqrt(3.1 + 3.1)),
                DiffractionDataPoint(5.03, 4.1 - 4.1, math.sqrt(4.1 + 4.1)),
                DiffractionDataPoint(5.04, 6.1 - 6.1, math.sqrt(6.1 + 6.1)),
                DiffractionDataPoint(5.05, 5.1 - 5.1, math.sqrt(5.1 + 5.1))]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    dp1 = make_dp()
    dp1 -= dp2
    assert is_all_equal(dp1, dp3)

    with pytest.raises(TypeError) as e_info:
        dp1 *= dp2
        dp1 /= dp2
        dp1 //= dp2
        dp1 += (1, 2, 3)
        dp1 -= (1, 2, 3)

    with pytest.raises(ValueError) as e_info:
        dp1 += [1, 2, 3]
        dp1 -= [1, 2, 3]
