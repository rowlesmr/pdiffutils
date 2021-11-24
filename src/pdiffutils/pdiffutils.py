# -*- coding: utf-8 -*-
"""
Created on Sun May  2 19:45:57 2021

@author: 184277J
"""
from __future__ import annotations

import math
import os
import copy
import operator
# import random
# from timeit import default_timer as timer  # use as start = timer() ...  end = timer()
from math import sqrt
from typing import List, Union, Tuple


def val_err_str(val: float, err: float) -> str:
    """
    Get a float representation of a value/error pair and create a string representation
    12.345 +/- 1.23 --> 12.3(12)
    12.345 +/- 0.012 -> 12.345(12
    12345 +/- 654  ---> 12300(650)
    :param val: float representing the value
    :param err: float representing the error in the value
    :return: a string representation of the value/error pair
    """
    err_sig_figs = 2  # future upgrade path is to allow user to set this
    dps = 2 - err_sig_figs
    if err < 10:
        while err < 10.:
            err *= 10
            dps += 1
        err = round(err, 0)
    else:  # err > 10
        while err > 100.:
            err /= 10
            dps -= 1
        err = round(err, 0) * 10 ** (-dps)
    val = round(val, dps)
    return f"{val:.{max(0, dps)}f}({err:.0f})"


class DiffractionDataPoint:
    """
    This is an X-ray data point.

    Each entry has a position, intensity and an error associated with the intensity.

    The position variable is named "angle", but this may also be energy (for ED) or time (for TOF).
    The name is just semantically related for "normal" diffraction
    """

    def __init__(self, x: float, y: float, error: float = None):
        self.x = x
        self.y = y
        self.e = math.sqrt(abs(self.y)) if error is None else abs(error)

    def negate(self):
        self.x *= -1.0

    def zeroOffset(self, offset: float):
        self.x += offset

    def average_with(self, ddps: List[DiffractionDataPoint], in_place: bool = True) -> DiffractionDataPoint:
        if not ddps:  # empty list
            return copy.deepcopy(self)

        ddp_mean = DiffractionDataPointMean()
        ddp_mean.increment(self)
        for ddp in ddps:
            if ddp.x != self.x:
                raise ValueError(f"trying to average data at {ddp.x} with data at {self.x}.")
            ddp_mean.increment(ddp)

        if in_place:
            self.x = ddp_mean.get_result().x
            self.y = ddp_mean.get_result().y
            self.e = ddp_mean.get_result().e
        else:
            return ddp_mean.get_result()

    def __str__(self) -> str:
        return f"{self.x:.5f} {self.y:.3f} {self.e:.3f}"

    def __repr__(self):
        return f"XRayDataPoint({self.x}, {self.y}, {self.e})"

    def __lt__(self, other: DiffractionDataPoint):
        if isinstance(other, DiffractionDataPoint):
            return self.x < other.x
        else:
            return NotImplemented

    def __le__(self, other: DiffractionDataPoint):
        if isinstance(other, DiffractionDataPoint):
            return self.x <= other.x
        else:
            return NotImplemented

    def __eq__(self, other: DiffractionDataPoint):
        if isinstance(other, DiffractionDataPoint):
            return math.isclose(self.x, other.x)
        else:
            return NotImplemented

    def __ne__(self, other: DiffractionDataPoint):
        if isinstance(other, DiffractionDataPoint):
            return not math.isclose(self.x, other.x)
        else:
            return NotImplemented

    def __ge__(self, other: DiffractionDataPoint):
        if isinstance(other, DiffractionDataPoint):
            return self.x >= other.x
        else:
            return NotImplemented

    def __gt__(self, other: DiffractionDataPoint):
        if isinstance(other, DiffractionDataPoint):
            return self.x > other.x
        else:
            return NotImplemented

    def __hash__(self):
        return hash((self.x, self.y, self.e))

    def __mul__(self, other: float) -> DiffractionDataPoint:
        """
        Overriding * operator to mean multipling intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Multiplication with type {type(other)} is undefined.")

        return DiffractionDataPoint(self.x,
                                    self.y * other,
                                    self.e * other)

    def __imul__(self, other: float) -> DiffractionDataPoint:
        """
        Overriding *= operator to mean multipling intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            return NotImplemented

        self.y *= other
        self.e *= other

        return self

    def __truediv__(self, other: float) -> DiffractionDataPoint:
        """
        Overriding / operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        return DiffractionDataPoint(self.x,
                                    self.y / other,
                                    self.e / other)

    def __itruediv__(self, other: float) -> DiffractionDataPoint:
        """
        Overriding /= operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            return NotImplemented

        self.y /= other
        self.e /= other
        return self

    def __floordiv__(self, other: float) -> DiffractionDataPoint:
        """
        Overriding // operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        return DiffractionDataPoint(self.x,
                                    self.y // other,
                                    self.e // other)

    def __ifloordiv__(self, other: float) -> DiffractionDataPoint:
        """
        Overriding //= operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            return NotImplemented

        self.y //= other
        self.e //= other
        return self

    def _add_sub(self, other: Union[float, DiffractionDataPoint], op: operator) -> DiffractionDataPoint:
        """
        helper method to override + and -

        Parameters
        ----------
        other : float, int, or DiffractionDataPoint
            DESCRIPTION.
        op : operator.add or operator.subtract
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if isinstance(other, (float, int)):
            return DiffractionDataPoint(self.x,
                                        op(self.y, other),
                                        self.e)

        if not isinstance(other, DiffractionDataPoint):
            raise TypeError(f"Addition/subtraction with type {type(other)} is undefined.")

        if self == other:
            return DiffractionDataPoint(self.x,
                                        op(self.y, other.y),
                                        math.sqrt(self.e ** 2 + other.e ** 2))
        else:
            raise ValueError(f"Angles are not equal: {self.x} & {other.x}.")

    def __add__(self, other: Union[float, DiffractionDataPoint]) -> DiffractionDataPoint:
        """
        Overriding + operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Parameters
        ----------
        other : DiffractionDataPoint.

        Returns
        -------
        DiffractionDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._add_sub(other, operator.add)

    def __sub__(self, other: Union[float, DiffractionDataPoint]) -> DiffractionDataPoint:
        """
        Overriding - operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        DiffractionDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._add_sub(other, operator.sub)

    def _iadd_isub(self, other: Union[float, DiffractionDataPoint], op: operator) -> DiffractionDataPoint:
        """
        helper function
        Overriding  += and -=  operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Ths does it in-place


        Parameters
        ----------
        other : DiffractionDataPoint, float or int
        op: operator.iadd or operator isub

        Returns
        -------
        DiffractionDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if isinstance(other, (float, int)):
            self.y = op(self.y, other)
            return self

        if not isinstance(other, DiffractionDataPoint):
            return NotImplemented

        if self.x != other.x:
            raise ValueError(f"Angles are not equal: {self.x} & {other.x}.")
        self.y = op(self.y, other.y)
        self.e = math.sqrt(self.e ** 2 + other.e ** 2)
        return self

    def __iadd__(self, other: Union[float, DiffractionDataPoint]) -> DiffractionDataPoint:
        """
        Overriding += operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Ths does it in-place


        Parameters
        ----------
        other : DiffractionDataPoint

        Returns
        -------
        DiffractionDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._iadd_isub(other, operator.iadd)

    def __isub__(self, other: Union[float, DiffractionDataPoint]) -> DiffractionDataPoint:
        """
        Overriding -= operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Ths does it in-place


        Parameters
        ----------
        other : DiffractionDataPoint

        Returns
        -------
        DiffractionDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._iadd_isub(other, operator.isub)


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------


class DiffractionPattern:
    """
    This is a Diffraction Pattern.
    It is essentially a list of DiffractionDataPoints.
    It takes either a filename or a list of DiffractionDataPoints
    """

    def __init__(self, diffpat: List[DiffractionDataPoint] = None, filename: str = None, meta=None):
        """
        Read in inital data in order to make everything.
        Angles must be strictly increasing

        Parameters
        ----------
        diffpat : a list of DiffractionDataPoints
        filename: a string representing the filename of the diffraction data
        meta: any metainformation about the diffraction pattern.

        Returns
        -------
        None.
        """
        self.filename = filename
        self.diffpat = Read.read(filename) if diffpat is None else diffpat
        self.meta = meta  # any information at all about this diffraction pattern. Can be anything of any type.

        self.xs = [xdp.x for xdp in self.diffpat]
        self.ys = [xdp.y for xdp in self.diffpat]
        self.es = [xdp.e for xdp in self.diffpat]

        step_sizes = [j - i for i, j in zip(self.xs[:-1], self.xs[1:])]
        self.ave_step_size = sum(step_sizes) / len(step_sizes) if len(step_sizes) > 1 else 0

    def negate(self):
        for d in self.diffpat:
            d.negate()

    def reverse(self):
        self.diffpat.reverse()

    def zeroOffset(self, offset: Union[float, List[float]]):
        if isinstance(offset, float):
            for d in self.diffpat:
                d.zeroOffset(offset)
        elif isinstance(offset, list):
            if len(offset) != len(self.diffpat):
                raise ValueError(f"Need {len(self.diffpat)} values. Was given {len(offset)} instead.")
            for d, o in zip(self.diffpat, offset):
                d.zeroOffset(o)

    def getMinAngle(self):
        return self.diffpat[0].x

    def getMaxAngle(self):
        return self.diffpat[-1].x

    def getData(self):
        return copy.deepcopy(self.diffpat)

    def __len__(self) -> int:
        return len(self.diffpat)

    def __str__(self):
        return "".join(str(d) + "\n" for d in self.diffpat)

    def __repr__(self):
        s = "DiffractionPattern([\n"
        for d in self.diffpat:
            s += repr(d) + ",\n"
        s = f"{s[:-2]}\n], {self.filename}, {self.meta})"
        return s

    def trim(self, min_x: float = -180, max_x: float = 180, in_place: bool = True) -> None | DiffractionPattern:
        """
        Trims a diffraction pattern such that there exist no angles less
        than min_angle, and no angles greater than max_angle.

        Parameters
        ----------
        min_x : float, optional
            the minimum angle you want to see in your diffraction pattern. The default is -180.
        max_x : float, optional
            the maximum angle you want to see in your diffraction pattern. The default is 180.
        in_place: bool, alter self, or return a new Diffractionpattern?
        """
        dp = [xdp for xdp in self.diffpat if min_x <= xdp.x <= max_x]

        if in_place:
            self.diffpat = dp
        else:
            return DiffractionPattern(diffpat=dp)

    def sort(self, in_place: bool = True):
        """
        sort the diffraction pattern, based on angles
        """
        if in_place:
            self.diffpat.sort()
        else:
            dp = copy.deepcopy(self.diffpat)
            dp.sort()
            return DiffractionPattern(diffpat=dp)

    def downsample(self, ds: int, in_place: bool = True) -> None | DiffractionPattern:
        """
        Downsamples the number of angles by averaging them.
        ds == 2 gives half the number of angles, 3 gives one third, and so on.
        Parameters
        ----------
        ds an integer describing the factor by which to downsample.

        Returns
        -------
        new downsampled diffraction pattern
        """
        a = self.xs
        i = self.ys
        e = self.es
        r = []

        for j in range(0, len(a), ds):
            a_new = sum(a[j:j + ds]) / ds
            i_new = sum(i[j:j + ds]) / ds
            e_new = sqrt(sum(map(lambda k: k * k, e[j:j + ds]))) / ds
            r.append(DiffractionDataPoint(a_new, i_new, e_new))

        if in_place:
            self.diffpat = r
        else:
            return DiffractionPattern(r)

    def split_on_zero(self) -> tuple[DiffractionPattern | None, DiffractionPattern | None]:
        """
        Gets a diffraction pattern that contains -ve and +ve angles, and returns two diffraction patterns:
        one from the +ve bit, and a negated version from the -ve side
        Returns
        -------
        A tuple containing two diffraction patterns. The first is from the +ve side, the second from the -ve side
        """
        min_ = self.getMinAngle()
        max_ = self.getMaxAngle()
        do_pos = False  # which angles do I do? positive? negative?
        do_neg = False
        pos = None
        neg = None

        if min_ < 0 and max_ > 0:
            do_pos = True
            do_neg = True
        elif min_ < 0 and max_ < 0:
            do_neg = True
        elif min_ > 0 and max_ > 0:
            do_pos = True

        if do_pos:
            pos = copy.deepcopy(self)
            pos.trim(0.0001, 180)
        if do_neg:
            neg = copy.deepcopy(self)
            neg.trim(-180, -0.0001)
            neg.negate()
            neg.sort()

        return pos, neg

    def average_with(self, dps: List[DiffractionPattern], in_place: bool = True) -> None | DiffractionPattern:
        # get deepcopy of the diffpats so I'm not plagued by pointer errors
        src = copy.deepcopy(self.diffpat)
        for dp in dps:
            src += copy.deepcopy(dp.diffpat)
        src.sort()

        dest = []
        tmp = []
        i = 0
        while i < len(src):
            p1 = src[i]
            j = i + 1
            while j < len(src):
                p2 = src[j]
                if p2 > p1:
                    break
                elif math.isclose(p2.x, p1.x):
                    tmp.append(p2)
                    j += 1
            dest.append(p1.average_with(tmp, in_place=False))
            tmp = []
            i = j

        if in_place:
            self.diffpat = dest
        else:
            return DiffractionPattern(diffpat=dest)

    def _add_and(self, other: DiffractionPattern, op: operator) -> DiffractionPattern:
        """
        This function is used by __add__ and __and__ to allow for + and &
        to mean summing or averaging diffraction patterns.

        The difference between the two functions is only a single operator, so
        to make it easier to maintain, I use the "operator" module to allow me
        to pass in a function that acts as an operator


        Parameters
        ----------
        other : DiffractionPattern.
        op: a function that acts as an operator

        Returns
        -------
        DiffractionPattern
            with the intensities summed or averaged together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.
        TypeError: if you try to use the wrong other type

        """
        if not isinstance(other, DiffractionPattern):
            raise TypeError(f"Addition with type {type(other)} is undefined.")

        # get deepcopy of the diffpats so I'm not plagued by pointer errors
        lst = copy.deepcopy(self.diffpat)
        lst_other = copy.deepcopy(other.diffpat)

        # concatenate the two diffpat list and sort in increasing angle
        lst = lst + lst_other
        lst.sort()

        # the list of XRayDataPoints is now sorted in ascending order
        # an XDP == XDP iff the angles are equal
        i = 0
        j = 0
        r = []

        # This gets a point (p1) in the sorted, concatenated list and compares it with the next point (p2)
        #   if the angle of p2 is greater than p1, we append p1 to the return list and move to the next point
        #   if the angle of p2 is the same as p1, we add the two points together, and assign it back to p1.
        #     p2 is then removed from the list, and the next point is assigned to p2 and the check happens again
        # There is an "operator" function in there, as the process is the same for + and &, so it's easier to
        #  maintain one set of code, and use appropriate operators in __add__ and __and__
        while i < len(lst):
            p1 = lst[i]
            j = i + 1
            while j < len(lst):
                p2 = lst[j]
                if p2 > p1:
                    break
                elif p2 == p1:
                    p1 = op(p1, p2)  # this is the bit that adds the the two DiffractionDataPoints
                    lst.pop(j)
            r.append(p1)
            i += 1
        return DiffractionPattern(diffpat=r)

    def _iadd_iand(self, other: DiffractionPattern, op: operator):
        """
        This function is used by __iadd__ and __iand__ to allow for + and &
        to mean summing or averaging diffraction patterns.

        The difference between the two functions is only a single operator, so
        to make it easier to maintain, I use the "operator" module to allow me
        to pass in a function that acts as an operator

        Parameters
        ----------
        other : DiffractionPattern.
        op: a function that acts as an operator

        Returns
        -------
        DiffractionPattern
            with the intensities summed or averaged together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.
        TypeError: if you try to use the wrong other type

        """
        if not isinstance(other, DiffractionPattern):
            return NotImplemented

        # concatenate the two diffpat list and sort in increasing angle
        self.diffpat = self.diffpat + other.diffpat
        self.diffpat.sort()

        # the list of XRayDataPoints is now sorted in ascending order
        # an XDP == XDP iff the angles are equal

        # This gets a point (p1) in the sorted, concatenated list and compares it with the next point (p2)
        #   if the angle of p2 is greater than p1, we append p1 to the return list and move to the next point
        #   if the angle of p2 is the same as p1, we add the two points together, and assign it back to p1.
        #     p2 is then removed from the list, and the next point is assigned to p2 and the check happens again
        # There is an "operator" function in there, as the process is the same for + and &, so it's easier to
        #  maintain one set of code, and use appropriate operators in __add__ and __and__
        i = 0
        j = 0
        while i < len(self.diffpat):
            p1 = self.diffpat[i]
            j = i + 1
            while j < len(self.diffpat):
                p2 = self.diffpat[j]
                if p2 > p1:
                    break
                elif p2 == p1:
                    p1 = op(p1, p2)  # this is the bit that adds the the two DiffractionDataPoints
                    self.diffpat.pop(j)
            i += 1
        return self

    def __mul__(self, other: float) -> DiffractionPattern:
        """
        Overriding * operator to mean multipling intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionPattern
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Multiplication with type {type(other)} is undefined.")

        r = copy.deepcopy(self.diffpat)

        for i in range(len(r)):
            r[i] = r[i] * other
        return DiffractionPattern(diffpat=r)

    def __truediv__(self, other: float) -> DiffractionPattern:
        """
        Overriding / operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionPattern
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        r = copy.deepcopy(self.diffpat)

        for i in range(len(r)):
            r[i] = r[i] / other
        return DiffractionPattern(diffpat=r)

    def __floordiv__(self, other: float) -> DiffractionPattern:
        """
        Overriding // operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionPattern
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        r = copy.deepcopy(self.diffpat)

        for i in range(len(r)):
            r[i] = r[i] // other
        return DiffractionPattern(diffpat=r)

    def __add__(self, other: Union[float, DiffractionPattern]) -> DiffractionPattern:
        """
        Overriding + operator to mean adding diffraction pattern intensities together if the same angle
        or adding a float or int to all intensities

        Parameters
        ----------
        other : DiffractionPattern or int or float

        Returns
        -------
        DiffractionPattern
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if isinstance(other, (float, int)):
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] + other
            return DiffractionPattern(diffpat=r)
        elif isinstance(other, list):
            if len(other) != len(self.diffpat):
                raise ValueError(
                    f"You gave a list of {len(other)} to subtract from {len(self)} intensities. They need to be the same.")
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] + other[i]
            return DiffractionPattern(diffpat=r)
        else:
            return self._add_and(other, operator.add)

    def __sub__(self, other: Union[float, DiffractionPattern]):
        """
        Overriding - operator to mean adding diffraction pattern intensities together if the same angle
        or adding a float or int to all intensities

        Parameters
        ----------
        other : DiffractionPattern or int or float

        Returns
        -------
        DiffractionPattern
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if isinstance(other, (float, int)):
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] - other
            return DiffractionPattern(diffpat=r)
        elif isinstance(other, list):
            if len(other) != len(self.diffpat):
                raise ValueError(
                    f"You gave a list of {len(other)} to subtract from {len(self.diffpat)} intensities. They need to be the same.")
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] - other[i]
            return DiffractionPattern(diffpat=r)
        else:
            return self._add_and(other, operator.sub)

    def __iadd__(self, other):
        """
        Overriding += operator to mean adding diffraction pattern intensities together if the same angle
        or adding a float or int to all intensities

        doing it in place

        Parameters
        ----------
        other : DiffractionPattern or int or float

        Returns
        -------
        DiffractionPattern
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if isinstance(other, (float, int)):
            for i in range(len(self.diffpat)):
                self.diffpat[i] += other
            return self
        elif isinstance(other, list):
            if len(other) != len(self.diffpat):
                raise ValueError(
                    f"You gave a list of {len(other)} to add to {len(self.diffpat)} intensities. They need to be the same.")
            for i in range(len(self.diffpat)):
                self.diffpat[i] += other[i]
            return self
        else:
            return self._iadd_iand(other, operator.iadd)

    def __isub__(self, other):
        """
        Overriding -= operator to mean adding diffraction pattern intensities together if the same angle
        or adding a float or int to all intensities or a list of floats or ints

        doing it in place

        Parameters
        ----------
        other : DiffractionPattern or int or float

        Returns
        -------
        DiffractionPattern
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if isinstance(other, (float, int)):
            for i in range(len(self.diffpat)):
                self.diffpat[i] -= other
            return self
        elif isinstance(other, list):
            if len(other) != len(self.diffpat):
                raise ValueError(
                    f"You gave a list of {len(other)} to subtract from {len(self.diffpat)} intensities. They need to be the same.")
            for i in range(len(self.diffpat)):
                self.diffpat[i] -= other[i]
            return self
        else:
            return self._iadd_iand(other, operator.isub)

    def __itruediv__(self, other):
        """
        Overriding /= operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionPattern
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        for i in range(len(self.diffpat)):
            self.diffpat[i] /= other
        return self

    def __ifloordiv__(self, other):
        """
        Overriding //= operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        DiffractionPattern
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, (float, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        for i in range(len(self.diffpat)):
            self.diffpat[i] //= other
        return self


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class InterpolatedDiffractionPattern(DiffractionPattern):

    def __init__(self, diffpat: List[DiffractionDataPoint], interp_stepsize: float, filename: str = None, meta=None):

        """
        Read in inital data in order to make everything.

        the angular range of the interpolation is determined by the data itself
          I try to go as wide as possible If you want to limit the angles, you can trim
          them later.

        Parameters
        ----------
        data : a string representing a filename, or a list of XRayDataPoints
        interp_stepsize: a float giving the stepsize I want in the inerpolation

        Returns
        -------
        None.

        """
        super().__init__(diffpat, filename, meta)
        self.diffpat = InterpolatedDiffractionPattern.spline(interp_stepsize, self)

    @staticmethod
    def spline(step_spline, dp: DiffractionPattern):
        # print(f"Now splining {dp.filename} in steps of {step_spline}")
        IDP = InterpolatedDiffractionPattern
        # these will hold the full splined data once I've finished
        x_spline = []
        y_spline = []
        e_spline = []

        # These are the original data that I am going to spline
        xs = dp.xs
        ys = dp.ys
        es = dp.es

        # arbitrary choice of cut-off size -  If i encounter a stepsize
        #   greater then 5*the average, then I'll assume that is the next module
        step_threshold = 5 * dp.ave_step_size

        # Now I start the spline process
        start_index = 0
        for stop_index, i in enumerate(range(1, len(xs)), start=1):
            step = xs[i] - xs[i - 1]
            if step >= step_threshold or i == len(xs) - 1:
                # this is the range of data I want to interpolate
                #  I want to trim off a few datapoints either side of the module edge
                #  so I don't have to deal with their noisy edges.
                DROP_POINTS = 5
                x_list = xs[start_index + DROP_POINTS:stop_index - DROP_POINTS]
                y_list = ys[start_index + DROP_POINTS:stop_index - DROP_POINTS]
                e_list = es[start_index + DROP_POINTS:stop_index - DROP_POINTS]

                # this is the interpolated data
                x_interp = IDP.generateInterpList(x_list[0], x_list[-1], step_spline)
                y_interp = IDP.cubic_interp1d(x_interp, x_list, y_list)
                e_interp = IDP.cubic_interp1d(x_interp, x_list, e_list)

                # put the just-interpolated-data into the whole list of data to be used to make the new DP
                x_spline += x_interp
                y_spline += y_interp
                e_spline += e_interp

                start_index = stop_index

        return [DiffractionDataPoint(x_spline[i], y_spline[i], e_spline[i]) for i in range(len(x_spline))]

    @staticmethod
    def spline2(interp, dp):
        """
        Interpolates a given DiffractionPattern using the provided interp array.
        Intensities _and_ errors are interpolated. This second one really should
        be looked at properly to see how best to do it. But, for now, it works.

        Parameters
        ----------
        interp : float list
            what angles to I want to interpolate
        dp : DiffractionPattern
            the DiffractionPattern I want to interpolate

        Raises
        ------
        ValueError
            if the provided interpolation list contains angles
            outside of the range of the data.

        Returns
        -------
        r : list of XRayDataPoints
            The interpolated data

        """
        angle = dp.xs
        intensity = dp.ys
        error = dp.es

        # check limits on interpolation
        if interp[0] < angle[0]:  # ie I'm interpolating before the data starts
            raise ValueError(f"Trying to interpolate out of range: {interp[0]} < {angle[0]}")
        elif interp[-1] > angle[-1]:  # ie I'm interpolating after the data ends
            raise ValueError(f"Trying to interpolate out of range: {interp[-1]} > {angle[-1]}")
        # interpolate intensity and error
        intensity_interp = InterpolatedDiffractionPattern.cubic_interp1d(interp, angle, intensity)
        error_interp = InterpolatedDiffractionPattern.cubic_interp1d(interp, angle, error)

        return [
            DiffractionDataPoint(interp[i], intensity_interp[i], error_interp[i])
            for i in range(len(interp))
        ]

    @staticmethod
    def generateInterpList(min_start, max_stop, step):
        """
        Generate a list to use as the interpolation array.

        The steps are generated from 0, with m

        Step size is always positive
        min_start < max_stop

        Parameters
        ----------
        min_start : float
            The minimum value the start angle could take on
        step : float, positive
            the step size between angles.
        max_stop : float > min_start
            The maximum value the stop angle could take on

        Returns
        -------
        A list of floats from ~min_start to ~max_stop in steps of step.
        The actual start and stop angles will always be in the range
        [min_start, max_stop]

        Raises
        ------

        ValueError is max_stop < min_start

        """
        if max_stop <= min_start:
            raise ValueError("max_stop is less than min_start. This list must be strictly increasing.")

        angle = 0.0
        i = 0

        # find the starting angle
        if min_start >= 0:
            while angle < min_start:
                angle = i * step
                i += 1
            if max_stop < min_start:
                angle -= step
        else:
            while angle > min_start:
                angle = -i * step
                i += 1
            if max_stop > min_start:
                angle += step

        # build the list to return
        lst = []
        i = angle // step + 1
        if max_stop > min_start:
            while angle < max_stop:
                lst.append(angle)
                angle = i * step
                i += 1
            angle -= step
        else:
            while angle > max_stop:
                lst.append(angle)
                angle = -i * step
                i += 1
        return lst

    @staticmethod
    def cubic_interp1d(x0: Union[float, List[float]], x: List[float], y: List[float], do_checks: bool = True):
        """
        Interpolate a 1-D function using cubic splines.
          x0 : a float or 1d-list of floats to interpolate at
          x  : a 1-D list of floats sorted in increasing order
          y  : a 1-D list of floats. The length of y must be equal to the length of x.

        Implement a trick to generate at first step the cholesky matrice L of
        the tridiagonal matrice A (thus L is a bidiagonal matrice that
        can be solved in two distinct loops).

        additional ref: www.math.uh.edu/~jingqiu/math4364/spline.pdf
        # original function code at: https://stackoverflow.com/a/48085583/36061

        This function is licenced under: Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)
        https://creativecommons.org/licenses/by-sa/3.0/
        Original Author raphael valentin
        Date 3 Jan 2018

        Modifications made to remove numpy dependencies; all sub-functions by MR

        This function, and all sub-functions, are licenced under: Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)

        Mod author: Matthew Rowles
        Date 3 May 2021
        """
        if do_checks:
            if len(x) != len(y):
                raise ValueError(f"Length of x ({len(x)}) != length of y ({len(y)}).")

            for i in range(1, len(x)):
                if x[i - 1] >= x[i]:
                    raise ValueError("x is not in strictly increasing order.")

        def diff(lst: List[float]) -> List[float]:
            """
            numpy.diff with default settings
            """
            length = len(lst) - 1
            r = [0.0] * length
            for k in range(length):
                r[k] = lst[k + 1] - lst[k]
            return r

        def list_searchsorted(listToInsert: List[float], insertInto: List[float]) -> List[int]:
            """
            numpy.searchsorted with default settings
            """

            def float_searchsorted(floatToInsert: float, insertInto: List[float]) -> int:
                """
                Helper function
                """
                for i in range(len(insertInto)):
                    if floatToInsert <= insertInto[i]:
                        return i
                return len(insertInto)

            return [float_searchsorted(item, insertInto) for item in listToInsert]

        def clip(lst: List[float], min_val: float, max_val: float, in_place: bool = False) -> List[float]:
            """
            numpy.clip
            """
            if not in_place:
                lst = lst[:]
            for k in range(len(lst)):
                if lst[k] < min_val:
                    lst[k] = min_val
                elif lst[k] > max_val:
                    lst[k] = max_val
            return lst

        def subtract(a: float, b: float) -> float:
            """
            returns a - b
            """
            return a - b

        if type(x0) is float:
            x0 = [x0]

        dim = len(x)

        xdiff = diff(x)
        ydiff = diff(y)

        # allocate buffer matrices
        Li: list = [0] * dim
        Li_1: list = [0] * (dim - 1)
        z: list = [0] * dim

        # fill diagonals Li and Li-1 and solve [L][y] = [B]
        Li[0] = sqrt(2 * xdiff[0])
        Li_1[0] = 0.0
        B0 = 0.0  # natural boundary
        z[0] = B0 / Li[0]

        for i in range(1, dim - 1, 1):
            Li_1[i] = xdiff[i - 1] / Li[i - 1]
            Li[i] = sqrt(2 * (xdiff[i - 1] + xdiff[i]) - Li_1[i - 1] * Li_1[i - 1])
            Bi = 6 * (ydiff[i] / xdiff[i] - ydiff[i - 1] / xdiff[i - 1])
            z[i] = (Bi - Li_1[i - 1] * z[i - 1]) / Li[i]

        i = dim - 1
        Li_1[i - 1] = xdiff[-1] / Li[i - 1]
        Li[i] = sqrt(2 * xdiff[-1] - Li_1[i - 1] * Li_1[i - 1])
        Bi = 0.0  # natural boundary
        z[i] = (Bi - Li_1[i - 1] * z[i - 1]) / Li[i]

        # solve [L.T][x] = [y]
        i = dim - 1
        z[i] = z[i] / Li[i]
        for i in range(dim - 2, -1, -1):
            z[i] = (z[i] - Li_1[i - 1] * z[i + 1]) / Li[i]

        # find index
        index = list_searchsorted(x0, x)
        index = clip(index, 1, dim - 1)

        xi1 = [x[num] for num in index]
        xi0 = [x[num - 1] for num in index]
        yi1 = [y[num] for num in index]
        yi0 = [y[num - 1] for num in index]
        zi1 = [z[num] for num in index]
        zi0 = [z[num - 1] for num in index]
        hi1 = list(map(subtract, xi1, xi0))

        # calculate cubic - all element-wise multiplication
        f0 = [0] * len(hi1)
        for j in range(len(f0)):
            f0[j] = zi0[j] / (6 * hi1[j]) * (xi1[j] - x0[j]) ** 3 + \
                    zi1[j] / (6 * hi1[j]) * (x0[j] - xi0[j]) ** 3 + \
                    (yi1[j] / hi1[j] - zi1[j] * hi1[j] / 6) * (x0[j] - xi0[j]) + \
                    (yi0[j] / hi1[j] - zi0[j] * hi1[j] / 6) * (xi1[j] - x0[j])

        # let's me return a float if a float was put in as x0
        if len(f0) == 1:
            f0 = f0[0]
        return f0


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class DiffractionExperiment:
    """
    This is a Diffraction Experiment.
    It is essentially a list of DiffractionPatterns.
    """

    def __init__(self, diffpats: List[DiffractionPattern], filename: str = None, meta=None):
        """
        Read in inital data in order to make everything.
        Angles must be strictly increasing

        Parameters
        ----------
        data : a string representing a filename, or a list of XRayDataPoints

        Returns
        -------
        None.
        """
        self.filename = filename
        self.diffpats = diffpats
        self.meta = meta  # any information at all about this diffraction pattern. Can be anything of any type.

    def negate(self):
        for d in self.diffpats:
            d.negate()

    def reverse(self):
        for d in self.diffpats:
            d.reverse()

    def zeroOffset(self, offset: float):
        for d in self.diffpats:
            d.zeroOffset(offset)

    def getData(self):
        return copy.deepcopy(self.diffpats)

    def __len__(self):
        return len(self.diffpats)

    def __str__(self):
        return "".join(str(d) + "\n" for d in self.diffpats)

    def __repr__(self):
        s = "DiffractionExperiment([\n"
        for d in self.diffpats:
            s += repr(d) + ",\n"
        s = f"{s[:-2]}\n], {self.filename}, {self.meta})"
        return s

    def trim(self, min_angle: Union[float, List[float]] = -180,
             max_angle: Union[float, List[float]] = 180, in_place: bool = True) -> None | DiffractionExperiment:
        """
        Trims a diffraction pattern such that there exist no angles less
        than min_angle, and no angles greater than max_angle.

        Parameters
        ----------
        min_angle : float, optional
            the minimum angle you want to see in your diffraction pattern. The default is -180.
        max_angle : float, optional
            the maximum angle you want to see in your diffraction pattern. The default is 180.
        """
        dps = []
        if isinstance(min_angle, float) and isinstance(max_angle, float):
            dps = [dp.trim(min_angle, max_angle, in_place=True) for dp in self.diffpats]
        elif isinstance(min_angle, float) and isinstance(max_angle, list):
            dps = [dp.trim(min_angle, max_x, in_place=True) for dp, max_x in zip(self.diffpats, max_angle)]
        elif isinstance(min_angle, list) and isinstance(max_angle, float):
            dps = [dp.trim(min_x, max_angle, in_place=True) for dp, min_x in zip(self.diffpats, min_angle)]
        elif isinstance(min_angle, list) and isinstance(max_angle, list):
            dps = [dp.trim(min_x, max_x, in_place=True) for dp, min_x, max_x in zip(self.diffpats, min_angle, max_angle)]

        if in_place:
            self.diffpats = dps
        else:
            return DiffractionExperiment(diffpats=dps)

    def sort(self, in_place: bool = True) -> None | DiffractionExperiment:
        """
        In-place sort of the diffraction pattern, based on angles
        """
        dps = [dp.sort(in_place=in_place) for dp in self.diffpats]
        if in_place:
            self.diffpats = dps
        else:
            return DiffractionExperiment(diffpats=dps)

    def downsample(self, ds: int, in_place: bool = True) -> None | DiffractionExperiment:
        """
        Downsamples the number of angles by averaging them.
        ds == 2 gives half the number of angles, 3 gives one third, and so on.
        Parameters
        ----------
        ds an integer describing the factor by which to downsample.
        """
        dps = [dp.downsample(ds, in_place=in_place) for dp in self.diffpats]
        if in_place:
            self.diffpats = dps
        else:
            return DiffractionExperiment(diffpats=dps)

    def average_patterns(self, num: int, is_rolling: bool = True, in_place: bool = True) -> None | DiffractionExperiment:
        """
        Averages num patterns. If rolling average, then
        patterns 1..num are averaged, then 2..num+1, 3--num+2 and so on.
        If not rolling average, then 1..num, num+1..num+num and so on.
        :param num: number of patterns to average
        :param is_rolling: do a rolling average?
        :return:
        """
        range_step = 1 if is_rolling else num
        dps = []

        for i in range(0, len(self.diffpats) - num, range_step):
            dp = copy.deepcopy(self.diffpats[i])
            tmp = [self.diffpats[i + j] for j in range(1, num)]
            dps.append(dp.average_with(tmp))

        if in_place:
            self.diffpats = dps
        else:
            return DiffractionExperiment(diffpats=dps)

    def interpolate(self, step: float, new_x: List[float] = None):
        pass

    def split_on_zero(self):
        """
        Gets a diffraction pattern that contains -ve and +ve angles, and returns two diffraction patterns:
        one from the +ve bit, and a negated version from the -ve side
        Returns
        -------
        A tuple containing two diffraction patterns. The first is from the +ve side, the second from the -ve side
        """
        pass


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class DiffractionDataPointMean:
    """
    This is a way to average a whole bunch of diffraction points by adding them one-by-one to a list
    and progressively calculating the average. You don't need to know how many you've added; this
    keeps track of it. Copied from https://commons.apache.org/proper/commons-math/javadocs/api-2.2/src-html/org/apache/commons/math/stat/descriptive/moment/FirstMoment.html#line.76
    """

    def __init__(self):
        self.n: int = 0  # counter to keep track of how many items you want to average
        self.x: float = float("nan")  # the x value of the ddp
        self.m1_y: float = float("nan")  # the running average of the intensity
        self.dev_y: float = float("nan")  # the difference between the value you're adding and the current average
        self.nDev_y: float = float("nan")  # the differences divided by n
        self.m1_e: float = float("nan")  # the running average of the square of the error
        self.dev_e: float = float("nan")  # the difference between the value you're adding and the current average
        self.nDev_e: float = float("nan")  # the differences divided by n

    def clear(self):
        self.n = 0
        self.x = float("nan")
        self.m1_y = float("nan")
        self.dev_y = float("nan")
        self.nDev_y = float("nan")
        self.m1_e = float("nan")
        self.dev_e = float("nan")
        self.nDev_e = float("nan")

    def get_result(self) -> DiffractionDataPoint:
        return DiffractionDataPoint(self.x, self.m1_y, math.sqrt(self.m1_e / self.n))

    def increment(self, ddp: DiffractionDataPoint):
        if self.n == 0:
            self.x = ddp.x
            self.m1_y = 0.0
            self.m1_e = 0.0

        if self.x != ddp.x:
            raise ValueError(f"The given value of {ddp.x} is not equal to the first value, {self.x}.")

        self.n += 1

        self.dev_y = ddp.y - self.m1_y
        self.nDev_y = self.dev_y / self.n
        self.m1_y += self.nDev_y

        self.dev_e = ddp.e ** 2 - self.m1_e
        self.nDev_e = self.dev_e / self.n
        self.m1_e += self.nDev_e


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------


class Read:
    """
    Reads in data from file into a DiffractionPattern
    can read in XY, XYE, DAT files. Skips the line if it encounters a non-float thing.
    """

    @staticmethod
    def read(filename: str):
        """
        This guesses which read method to try based on filename.
        :param filename: string representing the contents of the file
        :return: a list of XRayDataPoints
        """
        if filename.endswith(".xye"):
            return Read.xye(filename)
        elif filename.endswith(".xy"):
            return Read.xy(filename)
        elif filename.endswith(".dat"):
            return Read.dat(filename)
        # getting here means that I don't know how to read the file
        raise ValueError(f"I don't know how to read {filename}.")

    @staticmethod
    def xye(filename: str, is_xy: bool = False, sort: bool = True):
        lst = []
        increasing = True
        # print(f"Now reading {filename}.")
        with open(filename) as f:
            print(f"Reading from {os.path.abspath(filename)}.")
            for line in f:
                s = line.split()  # splits on whitespace
                try:
                    angle = float(s[0])
                    intensity = float(s[1])
                    error = float(s[2]) if not is_xy else None
                except ValueError:
                    continue
                xdp = DiffractionDataPoint(angle, intensity, error)
                lst.append(xdp)

                # compare angles to see if strictly increasing
                if len(lst) >= 2 and increasing:
                    first_angle = lst[-2].x  # second-last value
                    second_angle = lst[-1].x  # last value
                    if first_angle >= second_angle:  # ie the angle went down, and not up
                        increasing = False
        if not increasing:
            if sort:
                lst.sort()
            else:
                raise ValueError("X-axis must be in monotonically increasing order")
        return lst

    @staticmethod
    def xy(filename: str):
        return Read.xye(filename, is_xy=True)

    @staticmethod
    def dat(filename: str):
        lst = []
        # print(f"Now reading {filename}.")
        with open(filename) as f:
            print(f"Reading from {os.path.abspath(filename)}.")
            # this first for is to skip the comments and grab the start step stop values
            for line in f:
                s = line.split()  # splits on whitespace
                try:
                    start = float(s[0])
                    step = float(s[1])
                    stop = float(s[2])
                    break
                except (ValueError, IndexError):
                    continue
            num_points = int(round((stop - start) / step, 5)) + 1

            # This second loop is to get the intensities. It doesn't matter how many are on each line.
            #   This just sucks them all up, line by line.
            i = 0
            for line in f:
                tokens = line.split()  # splits on whitespace
                for token in tokens:
                    try:
                        angle = start + i * step
                        intensity = float(token)
                        xdp = DiffractionDataPoint(angle, intensity)
                        lst.append(xdp)
                        i += 1
                    except ValueError:  # this will be triggered if token isn't a float, so the token just gets skipped
                        continue
        if i != num_points:
            raise ValueError(f"Received {i} points, should have {num_points}.")
        return lst


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class Write:
    """
    Writes data in the various data formats
    """

    @staticmethod
    def _get_padding(dp: DiffractionPattern):
        x_pad = int(math.log10(max(dp.xs))) + 2
        y_pad = int(math.log10(max(dp.ys))) + 2
        e_pad = max(int(math.log10(max(dp.es))), 0) + 2
        return x_pad, y_pad, e_pad

    @staticmethod
    def _get_formatted(d: DiffractionDataPoint, x_dp: int, y_dp: int, e_dp: int, x_pad: int, y_pad: int, e_pad: int):
        a = f"{d.x:{1 + x_pad + x_dp}.{x_dp}f}"
        i = f"{d.y:{1 + y_pad + y_dp}.{y_dp}f}"
        e = f"{d.e:{1 + e_pad + e_dp}.{e_dp}f}"
        return a, i, e

    @staticmethod
    def xye(dp: DiffractionPattern, filename: str,
            dp_angle: int = 5, dp_intensity: int = 3, dp_error: int = 3, is_xy: bool = False):
        """
        To write a nice XYE file to a file.

        dp_ refers to the number of decimal places you want to see for angle, intensity, or error
        lp_ refers to the left padding of the numbers, so that all the decimal points line up

        For example:
        angle    = 45.368754896
        dp_angle = 5
        lp_angle = 4 ie 3 spaces allocated for digits, and one for a negative sign (if it exists)
        format(angle,f"{1+4+5}.{5}f") = '  45.36875'

        :param dp:
        :param filename:
        :param dp_angle: How many decimal points do you want on the angle?. The default is 5.
        :param dp_intensity: How many decimal points do you want on the intensity?. The default is 3.
        :param dp_error: How many decimal points do you want on the error?. The default is 3.
        :param is_xy:
        :return:
        """
        print(f"Writing to {os.path.abspath(filename)}.")
        with open(filename, "w") as f:
            diffpat = dp.diffpat
            x_pad, y_pad, e_pad = Write._get_padding(dp)

            for d in diffpat:
                a, i, e = Write._get_formatted(d, dp_angle, dp_intensity, dp_error, x_pad, y_pad, e_pad)
                if is_xy:
                    f.write(f"{a}{i}\n")
                else:
                    f.write(f"{a}{i}{e}\n")

    @staticmethod
    def xy(dp: DiffractionPattern, filename: str, dp_angle: int = 5, dp_intensity: int = 3, dp_error: int = 3):
        Write.xye(dp, filename, dp_angle, dp_intensity, dp_error, is_xy=True)

    @staticmethod
    def dat(dp: DiffractionPattern, filename: str, dp_angle: int = 5, dp_intensity: int = 0):
        diffpat = dp.diffpat
        x_pad, y_pad, _ = Write._get_padding(dp)
        x_start = diffpat[0].x
        x_stop = diffpat[-1].x
        num_points = len(diffpat)
        x_step = (x_stop - x_start) / (num_points - 1)

        print(f"Writing to {os.path.abspath(filename)}.")
        with open(filename, "w") as f:
            f.write("\n")
            f.write("\n")
            f.write("\n")
            f.write("\n")

            x_start = f"{x_start:.{dp_angle}f}"
            x_step = f"{x_step:.{dp_angle}f}"
            x_stop = f"{x_stop:.{dp_angle}f}"
            f.write(f"{x_start} {x_step} {x_stop}\n")
            for k, d in enumerate(diffpat, start=1):
                _, i, _ = Write._get_formatted(d, dp_angle, dp_intensity, 3, x_pad, y_pad, 3)
                f.write(f"{i}")
                if k % 10 == 0:
                    f.write("\n")

    @staticmethod
    def cif(dp: DiffractionPattern, filename: str, dp_angle: int = 5,
            data_block: str = "pdiffutils", x_type: str = "_pd_meas_2theta_scan", y_type: str = "_pd_meas_intensity_total"):
        diffpat = dp.diffpat
        x_pad, _, _ = Write._get_padding(dp)
        ys = [val_err_str(y, e) for y, e in zip(dp.ys, dp.es)]
        max_len_y = 0
        for y in ys:
            max_len_y = max(max_len_y, len(y))

        print(f"Writing to {os.path.abspath(filename)}.")
        with open(filename, "w") as f:
            f.write(f"data_{data_block}\n")
            f.write("\tloop_\n")
            f.write(f"\t{x_type}\n")
            f.write(f"\t{y_type}\n")
            for d, y in zip(diffpat, ys):
                x, _, _ = Write._get_formatted(d, dp_angle, 3, 3, x_pad, 3, 3)
                f.write(f"\t{x} {y.rjust(max_len_y, ' ')}\n")


def main():
    dp1 = DiffractionPattern(filename=r"C:\Users\184277j\Documents\GitHub\pdiffutils\data\dp1.xy")
    dp2 = DiffractionPattern(filename=r"C:\Users\184277j\Documents\GitHub\pdiffutils\data\dp2.xy")
    dp3 = DiffractionPattern(filename=r"C:\Users\184277j\Documents\GitHub\pdiffutils\data\dp3.xy")

    # dp1                  dp2                   dp3
    # 5.00000 2.100 1.414  5.00000 3.200 1.732   5.00000 5.300 2.236
    # 5.01000 4.100 2.000  5.01000 4.200 2.000   5.01000 6.300 2.449
    # 5.02000 4.100 2.000  5.02000 5.200 2.236   5.02000 5.300 2.236
    # 5.03000 3.100 1.732  5.03000 5.200 2.236   5.03000 3.300 1.732
    # 5.04000 6.100 2.449  5.04000 7.200 2.646   5.04000 7.300 2.646

    # dp1 = DiffractionPattern(diffpat=[DiffractionDataPoint(5,2)])
    # dp2 = DiffractionPattern(diffpat=[DiffractionDataPoint(5,3)])
    # dp3 = DiffractionPattern(diffpat=[DiffractionDataPoint(5,5)])

    print(dp1.average_with([dp2, dp3]))


if __name__ == "__main__":
    main()
