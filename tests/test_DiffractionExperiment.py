import copy

import pytest

from pdiffutils import DiffractionExperiment, DiffractionPattern, DiffractionDataPoint
import math
import operator


# remember, I've overridden ==
def is_all_equal(d1: DiffractionExperiment, d2: DiffractionExperiment):  # pragma: no cover
    for dp1, dp2 in zip(d1.diffpats, d2.diffpats):
        for ddp1, ddp2 in zip(dp1.diffpat, dp2.diffpat):
            if not math.isclose(ddp1.x, ddp2.x):
                return False
            if not math.isclose(ddp1.y, ddp2.y):
                return False
            if not math.isclose(ddp1.e, ddp2.e):
                return False
    return True

def test_construction():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)

    de = DiffractionExperiment(diffpats=[dp1, dp2])

    assert de.diffpats == [dp1, dp2]


def test_negate():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de1 = DiffractionExperiment(diffpats=[dp1, dp2])
    de1.negate()

    diffpat3 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    dp3.negate()
    diffpat4 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp4 = DiffractionPattern(diffpat=diffpat4)
    dp4.negate()
    de2 = DiffractionExperiment(diffpats=[dp3, dp4])
    assert is_all_equal(de1, de2)


def test_reverse():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de1 = DiffractionExperiment(diffpats=[dp1, dp2])
    de1.reverse()

    diffpat3 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    dp3.reverse()
    diffpat4 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp4 = DiffractionPattern(diffpat=diffpat4)
    dp4.reverse()
    de2 = DiffractionExperiment(diffpats=[dp3, dp4])
    assert is_all_equal(de1, de2)


def test_zero_offset():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de1 = DiffractionExperiment(diffpats=[dp1, dp2])
    de1.zero_offset(0.1)

    diffpat3 = [DiffractionDataPoint(5.10, 2.1),
                DiffractionDataPoint(5.11, 4.1),
                DiffractionDataPoint(5.12, 4.1),
                DiffractionDataPoint(5.13, 3.1),
                DiffractionDataPoint(5.14, 6.1)]
    dp3 = DiffractionPattern(diffpat=diffpat3)

    diffpat4 = [DiffractionDataPoint(5.10, 4.2),
                DiffractionDataPoint(5.11, 2.2),
                DiffractionDataPoint(5.12, 4.2),
                DiffractionDataPoint(5.13, 3.2),
                DiffractionDataPoint(5.14, 5.2)]
    dp4 = DiffractionPattern(diffpat=diffpat4)
    de2 = DiffractionExperiment(diffpats=[dp3, dp4])
    assert is_all_equal(de1, de2)


def test_len():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de1 = DiffractionExperiment(diffpats=[dp1, dp2])

    assert len(de1) == 2


def test_trim():
    def make_de():
        diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                    DiffractionDataPoint(5.01, 4.1),
                    DiffractionDataPoint(5.02, 4.1),
                    DiffractionDataPoint(5.03, 3.1),
                    DiffractionDataPoint(5.04, 6.1),
                    DiffractionDataPoint(5.05, 7.1)]
        dp1 = DiffractionPattern(diffpat=diffpat1)
        diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                    DiffractionDataPoint(5.01, 2.2),
                    DiffractionDataPoint(5.02, 4.2),
                    DiffractionDataPoint(5.03, 3.2),
                    DiffractionDataPoint(5.04, 5.2),
                    DiffractionDataPoint(5.05, 6.2)]
        dp2 = DiffractionPattern(diffpat=diffpat2)
        return DiffractionExperiment(diffpats=[dp1, dp2])

    diffpat1 = [DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de2 = DiffractionExperiment(diffpats=[dp1,dp2])

    de1 = make_de()
    de = de1.trim(5.005, 5.045, in_place= False)
    assert is_all_equal(de, de2)
    de1.trim(5.005, 5.045, in_place=True)
    assert is_all_equal(de1, de2)

    diffpat1 = [DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de2 = DiffractionExperiment(diffpats=[dp1,dp2])

    de1 = make_de()
    de = de1.trim([5.015, 5.005], 5.045, in_place= False)
    assert is_all_equal(de, de2)
    de1.trim([5.015, 5.005], 5.045, in_place=True)
    assert is_all_equal(de1, de2)

    diffpat1 = [DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de2 = DiffractionExperiment(diffpats=[dp1, dp2])

    de1 = make_de()
    de = de1.trim(5.005, [5.035, 5.045], in_place=False)
    assert is_all_equal(de, de2)
    de1.trim(5.005, [5.035, 5.045], in_place=True)
    assert is_all_equal(de1, de2)

    diffpat1 = [DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de2 = DiffractionExperiment(diffpats=[dp1, dp2])

    de1 = make_de()
    de = de1.trim([5.005, 5.015], [5.035, 5.045], in_place=False)
    assert is_all_equal(de, de2)
    de1.trim([5.005, 5.015], [5.035, 5.045], in_place=True)
    assert is_all_equal(de1, de2)


def test_sort():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de1 = DiffractionExperiment(diffpats=[dp1, dp2])

    diffpat3 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1)]
    dp3 = DiffractionPattern(diffpat=diffpat3)

    diffpat4 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2)]
    dp4 = DiffractionPattern(diffpat=diffpat4)
    de2 = DiffractionExperiment(diffpats=[dp3, dp4])

    de = de1.sort(in_place=False)
    assert is_all_equal(de, de2)
    de1.sort(in_place=True)
    assert is_all_equal(de1, de2)


def test_downsample():
    diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                DiffractionDataPoint(5.01, 4.1),
                DiffractionDataPoint(5.02, 4.1),
                DiffractionDataPoint(5.03, 3.1),
                DiffractionDataPoint(5.04, 6.1),
                DiffractionDataPoint(5.05, 7.1)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                DiffractionDataPoint(5.01, 2.2),
                DiffractionDataPoint(5.02, 4.2),
                DiffractionDataPoint(5.03, 3.2),
                DiffractionDataPoint(5.04, 5.2),
                DiffractionDataPoint(5.05, 6.2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de1 = DiffractionExperiment(diffpats=[dp1, dp2])

    diffpat3 = [DiffractionDataPoint(5.01, (2.1 + 4.1 + 4.1) / 3, math.sqrt(2.1 + 4.1 + 4.1) / 3),
                DiffractionDataPoint(5.04, (3.1 + 6.1 + 7.1) / 3, math.sqrt(3.1 + 6.1 + 7.1) / 3)]
    dp3 = DiffractionPattern(diffpat=diffpat3)

    diffpat4 = [DiffractionDataPoint(5.01, (4.2 + 2.2 + 4.2) / 3, math.sqrt(4.2 + 2.2 + 4.2) / 3),
                DiffractionDataPoint(5.04, (3.2 + 5.2 + 6.2) / 3, math.sqrt(3.2 + 5.2 + 6.2) / 3)]
    dp4 = DiffractionPattern(diffpat=diffpat4)
    de2 = DiffractionExperiment(diffpats=[dp3, dp4])

    de = de1.downsample(3, in_place=False)
    assert is_all_equal(de, de2)
    de1.downsample(3, in_place=True)
    assert is_all_equal(de1, de2)


def test_average_patterns():
    def make_de():
        diffpat1 = [DiffractionDataPoint(5.00, 2.1),
                    DiffractionDataPoint(5.01, 4.1),
                    DiffractionDataPoint(5.02, 4.1),
                    DiffractionDataPoint(5.03, 3.1),
                    DiffractionDataPoint(5.04, 6.1),
                    DiffractionDataPoint(5.05, 7.1)]
        dp1 = DiffractionPattern(diffpat=diffpat1)
        diffpat2 = [DiffractionDataPoint(5.00, 4.2),
                    DiffractionDataPoint(5.01, 2.2),
                    DiffractionDataPoint(5.02, 4.2),
                    DiffractionDataPoint(5.03, 3.2),
                    DiffractionDataPoint(5.04, 5.2),
                    DiffractionDataPoint(5.05, 6.2)]
        dp2 = DiffractionPattern(diffpat=diffpat2)
        diffpat3 = [DiffractionDataPoint(5.00, 5.3),
                    DiffractionDataPoint(5.01, 3.3),
                    DiffractionDataPoint(5.02, 7.3),
                    DiffractionDataPoint(5.03, 3.3),
                    DiffractionDataPoint(5.04, 4.3),
                    DiffractionDataPoint(5.05, 2.3)]
        dp3 = DiffractionPattern(diffpat=diffpat3)
        diffpat4 = [DiffractionDataPoint(5.00, 3.4),
                    DiffractionDataPoint(5.01, 7.4),
                    DiffractionDataPoint(5.02, 8.4),
                    DiffractionDataPoint(5.03, 4.4),
                    DiffractionDataPoint(5.04, 6.4),
                    DiffractionDataPoint(5.05, 3.4)]
        dp4 = DiffractionPattern(diffpat=diffpat4)
        return DiffractionExperiment(diffpats=[dp1, dp2, dp3, dp4])

    #test rolling
    diffpat1 = [DiffractionDataPoint(5.00, (2.1+4.2)/2, math.sqrt(2.1+4.2)/2),
                DiffractionDataPoint(5.01, (4.1+2.2)/2, math.sqrt(4.1+2.2)/2),
                DiffractionDataPoint(5.02, (4.1+4.2)/2, math.sqrt(4.1+4.2)/2),
                DiffractionDataPoint(5.03, (3.1+3.2)/2, math.sqrt(3.1+3.2)/2),
                DiffractionDataPoint(5.04, (6.1+5.2)/2, math.sqrt(6.1+5.2)/2),
                DiffractionDataPoint(5.05, (7.1+6.2)/2, math.sqrt(7.1+6.2)/2)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, (4.2+5.3)/2, math.sqrt(4.2+5.3)/2),
                DiffractionDataPoint(5.01, (2.2+3.3)/2, math.sqrt(2.2+3.3)/2),
                DiffractionDataPoint(5.02, (4.2+7.3)/2, math.sqrt(4.2+7.3)/2),
                DiffractionDataPoint(5.03, (3.2+3.3)/2, math.sqrt(3.2+3.3)/2),
                DiffractionDataPoint(5.04, (5.2+4.3)/2, math.sqrt(5.2+4.3)/2),
                DiffractionDataPoint(5.05, (6.2+2.3)/2, math.sqrt(6.2+2.3)/2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    diffpat3 = [DiffractionDataPoint(5.00, (5.3+3.4)/2, math.sqrt(5.3+3.4)/2),
                DiffractionDataPoint(5.01, (3.3+7.4)/2, math.sqrt(3.3+7.4)/2),
                DiffractionDataPoint(5.02, (7.3+8.4)/2, math.sqrt(7.3+8.4)/2),
                DiffractionDataPoint(5.03, (3.3+4.4)/2, math.sqrt(3.3+4.4)/2),
                DiffractionDataPoint(5.04, (4.3+6.4)/2, math.sqrt(4.3+6.4)/2),
                DiffractionDataPoint(5.05, (2.3+3.4)/2, math.sqrt(2.3+3.4)/2)]
    dp3 = DiffractionPattern(diffpat=diffpat3)
    de2 = DiffractionExperiment(diffpats=[dp1, dp2, dp3])

    de1 = make_de()
    de = de1.average_patterns(2, is_rolling=True, in_place=False)
    assert is_all_equal(de, de2)
    de1.average_patterns(2, is_rolling=True, in_place=True)
    assert is_all_equal(de1, de2)

    #test not rolling
    diffpat1 = [DiffractionDataPoint(5.00, (2.1+4.2)/2, math.sqrt(2.1+4.2)/2),
                DiffractionDataPoint(5.01, (4.1+2.2)/2, math.sqrt(4.1+2.2)/2),
                DiffractionDataPoint(5.02, (4.1+4.2)/2, math.sqrt(4.1+4.2)/2),
                DiffractionDataPoint(5.03, (3.1+3.2)/2, math.sqrt(3.1+3.2)/2),
                DiffractionDataPoint(5.04, (6.1+5.2)/2, math.sqrt(6.1+5.2)/2),
                DiffractionDataPoint(5.05, (7.1+6.2)/2, math.sqrt(7.1+6.2)/2)]
    dp1 = DiffractionPattern(diffpat=diffpat1)
    diffpat2 = [DiffractionDataPoint(5.00, (5.3+3.4)/2, math.sqrt(5.3+3.4)/2),
                DiffractionDataPoint(5.01, (3.3+7.4)/2, math.sqrt(3.3+7.4)/2),
                DiffractionDataPoint(5.02, (7.3+8.4)/2, math.sqrt(7.3+8.4)/2),
                DiffractionDataPoint(5.03, (3.3+4.4)/2, math.sqrt(3.3+4.4)/2),
                DiffractionDataPoint(5.04, (4.3+6.4)/2, math.sqrt(4.3+6.4)/2),
                DiffractionDataPoint(5.05, (2.3+3.4)/2, math.sqrt(2.3+3.4)/2)]
    dp2 = DiffractionPattern(diffpat=diffpat2)
    de2 = DiffractionExperiment(diffpats=[dp1, dp2])

    de1 = make_de()
    de = de1.average_patterns(2, is_rolling=False, in_place=False)
    assert is_all_equal(de, de2)
    de1.average_patterns(2, is_rolling=False, in_place=True)
    assert is_all_equal(de1, de2)

@pytest.mark.skip
def test_interpolate():
    pass


@pytest.mark.skip
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
