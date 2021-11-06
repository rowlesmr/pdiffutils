# -*- coding: utf-8 -*-
"""
Created on Sun May  2 19:45:57 2021

@author: 184277J
"""

import math
import os
import copy
import operator
from math import sqrt


class XRayDataPoint:
    """
    This is an X-ray data point.

    Each entry has a position, intensity and an error associated with the intensity.

    The position variable is named "angle", but this may also be energy (for ED) or time (for TOF).
    The name is just semantically related for "normal" diffraction
    """

    def __init__(self, angle, intensity, error=None):
        self.angle = angle
        self.intensity = intensity
        if error is None:
            self.error = math.sqrt(abs(self.intensity))
        else:
            self.error = abs(error)

    def negate(self):
        self.angle *= -1.0

    def zeroOffset(self, offset):
        self.angle += offset

    def format_output(self, dp_angle=5, dp_intensity=3, dp_error=3, lp_angle=4, lp_intensity=6, lp_error=4):
        """
        Format the output of an XRayDataPoint to make it look nice.
        This is used primarily when printing DiffractionPatterns to a file.
        
        dp_ refers to the number of decimal places you want to see for angle, intensity, or error
        lp_ refers to the left padding of the numbers, so that all the decimal points line up
        
        For example:
        angle    = 45.368754896
        dp_angle = 5
        lp_angle = 4 ie 3 spaces allocated for digits, and one for a negative sign (if it exists)
        format(angle,f"{1+4+5}.{5}f") = '  45.36875'
        
        
        Parameters
        ----------
        dp_angle : int, optional
            How many decimal points do you want on the angle?. The default is 5.
        dp_intensity : int, optional
            How many decimal points do you want on the intensity?. The default is 3.
        dp_error : int, optional
            How many decimal points do you want on the error?. The default is 3.
        lp_angle : int, optional
            How many digits do you want to show on the left of the decimal point on the angle?. The default is 4.
        lp_intensity : int, optional
            How many digits do you want to show on the left of the decimal point on the intensity?. The default is 6.
        lp_error : int, optional
            How many digits do you want to show on the left of the decimal point on the error?. The default is 4.

        Returns
        -------
        str
            A nicely formatted string representing the values of the XRayDataPoint.

        """
        a = format(self.angle, f"{1 + lp_angle + dp_angle}.{dp_angle}f")
        i = format(self.intensity, f"{1 + lp_intensity + dp_intensity}.{dp_intensity}f")
        e = format(self.error, f"{1 + lp_error + dp_error}.{dp_error}f")
        return f"{a} {i} {e}"

    def __str__(self):
        return self.format_output(dp_angle=5, dp_intensity=3, dp_error=3)

    def __repr__(self):
        return f"XRayDataPoint({self.angle}, {self.intensity}, {self.error})"

    def __lt__(self, other):
        if isinstance(other, XRayDataPoint):
            return self.angle < other.angle
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, XRayDataPoint):
            return self.angle <= other.angle
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, XRayDataPoint):
            return math.isclose(self.angle, other.angle)
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, XRayDataPoint):
            return not math.isclose(self.angle, other.angle)
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, XRayDataPoint):
            return self.angle >= other.angle
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, XRayDataPoint):
            return self.angle > other.angle
        else:
            return NotImplemented

    def __hash__(self):
        return hash((self.angle, self.intensity, self.error))

    def __mul__(self, other):
        """
        Overriding * operator to mean multipling intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        XRayDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Multiplication with type {type(other)} is undefined.")

        return XRayDataPoint(self.angle,
                             self.intensity * other,
                             self.error * other)

    def __imul__(self, other):
        """
        Overriding *= operator to mean multipling intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        XRayDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not (isinstance(other, float) or isinstance(other, int)):
            return NotImplemented

        self.intensity *= other
        self.error *= other

        return self

    def __truediv__(self, other):
        """
        Overriding / operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        XRayDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        return XRayDataPoint(self.angle,
                             self.intensity / other,
                             self.error / other)

    def __itruediv__(self, other):
        """
        Overriding /= operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        XRayDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not (isinstance(other, float) or isinstance(other, int)):
            return NotImplemented

        self.intensity /= other
        self.error /= other
        return self

    def __floordiv__(self, other):
        """
        Overriding // operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        XRayDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        return XRayDataPoint(self.angle,
                             self.intensity // other,
                             self.error // other)

    def __ifloordiv__(self, other):
        """
        Overriding //= operator to mean dividing intensities by an int or float

        Parameters
        ----------
        other : int or float.

        Returns
        -------
        XRayDataPoint
            with the intensities scaled by other together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not (isinstance(other, float) or isinstance(other, int)):
            return NotImplemented

        self.intensity //= other
        self.error //= other
        return self

    def _add_sub(self, other, operator):
        """
        helper method to override + and -

        Parameters
        ----------
        other : float, int, or XRayDataPoint
            DESCRIPTION.
        operator : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if isinstance(other, float) or isinstance(other, int):
            return XRayDataPoint(self.angle,
                                 operator(self.intensity, other),
                                 self.error)

        if not isinstance(other, XRayDataPoint):
            raise TypeError(f"Addition/subtraction with type {type(other)} is undefined.")

        if self == other:
            return XRayDataPoint(self.angle,
                                 operator(self.intensity, other.intensity),
                                 math.sqrt(self.error ** 2 + other.error ** 2))
        else:
            raise ValueError(f"Angles are not equal: {self.angle} & {other.angle}.")

    def __add__(self, other):
        """
        Overriding + operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._add_sub(other, operator.add)

    def __sub__(self, other):
        """
        Overriding - operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._add_sub(other, operator.sub)

    def _iadd_isub(self, other, operator):
        """
        helper function 
        Overriding  += and -=  operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Ths does it in-place


        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if isinstance(other, float) or isinstance(other, int):
            self.intensity = operator(self.intensity, other)
            return self

        if not isinstance(other, XRayDataPoint):
            return NotImplemented

        if self.angle == other.angle:
            self.intensity = operator(self.intensity, other.intensity)
            self.error = math.sqrt(self.error ** 2 + other.error ** 2)
            return self
        else:
            raise ValueError(f"Angles are not equal: {self.angle} & {other.angle}.")

    def __iadd__(self, other):
        """
        Overriding += operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Ths does it in-place


        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._iadd_isub(other, operator.iadd)

    def __isub__(self, other):
        """
        Overriding -= operator to mean adding intensities together if the same angle
        Can also add a float or int. This will increase the intensity, but not alter the error

        Ths does it in-place


        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._iadd_isub(other, operator.isub)

    def __and__(self, other):
        """
        Overriding & operator to mean averaging intensities together if the same angle

        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint:
            with the intensities averaged together, and errors done nicely.

        Raises
        ------
        ValueError: if wrong type or angles not equal used.

        """
        if not isinstance(other, XRayDataPoint):
            raise TypeError(f"Averaging with type {type(other)} is undefined.")

        if self == other:
            return XRayDataPoint(self.angle,
                                 (self.intensity + other.intensity) / 2,
                                 math.sqrt(self.error ** 2 + other.error ** 2) / 2)
        else:
            raise ValueError(f"Angles are not equal: {self.angle} & {other.angle}")

    def __iand__(self, other):
        """
        Overriding &= operator to mean averaging intensities together if the same angle

        This does it inplace

        Parameters
        ----------
        other : XRayDataPoint.

        Returns
        -------
        XRayDataPoint
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        if not isinstance(other, XRayDataPoint):
            return NotImplemented

        if self == other:
            self.intensity = (self.intensity + other.intensity) / 2
            self.error = math.sqrt(self.error ** 2 + other.error ** 2) / 2
            return self
        else:
            raise ValueError(f"Angles are not equal: {self.angle} & {other.angle}.")


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class DiffractionPattern:
    """
    This is a Diffraction Pattern.

    It is essentially a list of XRayDataPoints.

    It takes either a filename or a list of XRayDataPoints

    """

    def __init__(self, data, aveStepSize=None, do_deep_copy=True):
        """
        Read in inital data in order to make everything.

        Angles must be strictly increasing


        Parameters
        ----------
        data : a string representing a filename, or a list of XRayDataPoints
        aveStepSize: a float containing the average step size, if already known.

        Returns
        -------
        None.

        Raises
        ------
        ValueError if angles don't strictly increase

        """

        self.diffpat = []
        self.aveStepSize = aveStepSize  # this will get used when splining. Eventually. When I get around to doing it.
        self.filename = None

        # This is for nice printing of the results - see writeToFile
        self.angle_padding = 0
        self.intensity_padding = 0
        self.error_padding = 0

        # deals with XY and XYE files with no comments and no missing lines
        if type(data) == str:  # treat it as a filename
            self.diffpat = self._readXYorXYE(data, aveStepSize)

        elif type(data) == list and all(isinstance(x, XRayDataPoint) for x in data):
            if do_deep_copy:
                self.diffpat = copy.deepcopy(data)
            else:
                self.diffpat = data

            n = 0
            tot = 0
            for i in range(1, len(self.diffpat)):
                n += 1
                tot += self.diffpat[i].angle - self.diffpat[i - 1].angle

                self.angle_padding = max(self.angle_padding, len(str(round(self.diffpat[i].angle))))
                self.intensity_padding = max(self.intensity_padding, len(str(round(self.diffpat[i].intensity))))
                self.error_padding = max(self.error_padding, len(str(round(self.diffpat[i].error))))

            if self.aveStepSize is None:
                self.aveStepSize = tot / n

        else:
            raise TypeError(f"Don't know how to read that input ({type(data)}).")

    def _readXYorXYE(self, filename, aveStepSize=None):
        lst = []
        n = 0
        tot = 0

        # print(f"Now reading {filename}.")
        with open(filename) as f:
            print(f"Reading from {os.path.abspath(filename)}.")
            self.filename = filename
            for line in f:
                s = line.split()  # splits on whitespace
                angle = float(s[0])
                intensity = float(s[1])
                if len(s) == 3:
                    error = float(s[2])
                else:
                    error = None

                xdp = XRayDataPoint(angle, intensity, error)

                self.angle_padding = max(self.angle_padding, len(str(round(xdp.angle))))
                self.intensity_padding = max(self.intensity_padding, len(str(round(xdp.intensity))))
                self.error_padding = max(self.error_padding, len(str(round(xdp.error))))

                lst.append(xdp)

                # compare angles to see if strictly increasing
                if len(lst) >= 2:
                    first_angle = lst[-2].angle  # second-last value
                    second_angle = lst[-1].angle  # last value

                    if first_angle >= second_angle:  # ie the angle went down, and not up
                        raise ValueError("Angles are not strictly increasing.")

                    n += 1
                    tot += second_angle - first_angle  # to calculate the averagestepsize

        if aveStepSize is None:
            self.aveStepSize = tot / n

        return lst

    def negate(self):
        for d in self.diffpat:
            d.negate()

    def reverse(self):
        self.diffpat.reverse()

    def zeroOffset(self, offset):
        for d in self.diffpat:
            d.zeroOffset(offset)

    def getMinAngle(self):
        return self.diffpat[0].angle

    def getMaxAngle(self):
        return self.diffpat[-1].angle

    def getNumOfDataPoints(self):
        return len(self.diffpat)

    def getAverageStepSize(self):
        return self.aveStepSize

    def getData(self):
        return copy.deepcopy(self.diffpat)

    def getAngles(self):
        """
        Get all the angles in the DP as a list

        Returns
        -------
        list of floats.

        """
        r = []
        for x in self.diffpat:
            r.append(x.angle)
        return r

    def getIntensities(self):
        """
        Get all the intensities in the DP as a list

        Returns
        -------
        list of floats.

        """
        r = []
        for x in self.diffpat:
            r.append(x.intensity)
        return r

    def getErrors(self):
        """
        Get all the erros in the DP as a list

        Returns
        -------
        list of floats.

        """
        r = []
        for x in self.diffpat:
            r.append(x.error)
        return r

    def __len__(self):
        return len(self.diffpat)

    def __str__(self):
        s = ""
        for d in self.diffpat:
            s += str(d) + "\n"
        return s

    def __repr__(self):
        s = "DiffractionPattern[\n"
        for d in self.diffpat:
            s += repr(d) + ",\n"
        s = s[:-2] + "\n]"
        return s

    def setPaddingPrettyPrinting(self, ang_, int_, err_):
        self.angle_padding = max(self.angle_padding, len(str(round(ang_))))
        self.intensity_padding = max(self.intensity_padding, len(str(round(int_))))
        self.error_padding = max(self.error_padding, len(str(round(err_))))

    def writeToFile(self, filename, dp_angle=5, dp_intensity=3, dp_error=3):
        print(f"Writing to {os.path.abspath(filename)}.")
        f = open(filename, "w")
        for d in self.diffpat:
            f.write(" " + d.format_output(dp_angle, dp_intensity, dp_error, self.angle_padding, self.intensity_padding,
                                          self.error_padding) + "\n")
        f.close()

    def trim(self, min_angle=-180, max_angle=180):
        """
        Trims a diffraction pattern such that there exist no angles less
        than min_angle, and no angles greater than max_angle.

        Parameters
        ----------
        min_angle : float, optional
            the minimum angle you want to see in your diffraction pattern. The default is -180.
        max_angle : float, optional
            the maximum angle you want to see in your diffraction pattern. The default is 180.

        Returns
        -------
        DiffractionPattern.

        """
        self.diffpat = [xdp for xdp in self.diffpat if xdp.angle >= min_angle and xdp.angle <= max_angle]

    def sort(self):
        """
        In-place sort of the diffraction pattern, based on angles

        Returns
        -------
        None.

        """
        self.diffpat.sort()

    def downsample(self, ds):
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
        a = self.getAngles()
        i = self.getIntensities()
        e = self.getErrors()
        r = []

        for j in range(0,len(a),ds):
            a_new = sum(a[j:j+ds]) / ds
            i_new = sum(i[j:j+ds]) / ds
            e_new = sqrt(sum(map(lambda i: i * i, e[j:j+ds]))) / ds
            r.append(XRayDataPoint(a_new, i_new, e_new))

        return DiffractionPattern(r)

    def split_on_zero(self):
        """
        Gets a diffraction pattern that contains -ve and +ve angles, and returns two diffraction patterns:
        one from the +ve bit, and a negated version from the -ve side
        Returns
        -------
        A tuple containing two diffraction patterns. The first is from the +ve side, the second from the -ve side
        """
        min_ = self.getMinAngle()
        max_ = self.getMaxAngle()
        do_pos = False # which angles do I do? positive? negative?
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

        return (pos, neg)
        


    def _add_and(self, other, operator):
        """

        This function is used by __add__ and __and__ to allow for + and &
        to mean summing or averging diffraction patterns.

        The difference between the two functions is only a single operator, so
        to make it easier to maintain, I use the "operator" module to allow me
        to pass in a function that acts as an operator


        Parameters
        ----------
        other : DiffractionPattern.
        operator: a function that acts as an operator

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
                    p1 = operator(p1, p2)  # this is the bit that adds the the two XRayDataPoints
                    lst.pop(j)
            r.append(p1)
            i += 1
        return DiffractionPattern(r,
                                  do_deep_copy=False)  # the initialisation takes a deepcopy of any incoming XRayDataPoint list

    def _iadd_iand(self, other, operator):
        """

        This function is used by __iadd__ and __iand__ to allow for + and &
        to mean summing or averaging diffraction patterns.

        The difference between the two functions is only a single operator, so
        to make it easier to maintain, I use the "operator" module to allow me
        to pass in a function that acts as an operator


        Parameters
        ----------
        other : DiffractionPattern.
        operator: a function that acts as an operator

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
                    p1 = operator(p1, p2)  # this is the bit that adds the the two XRayDataPoints
                    self.diffpat.pop(j)
            i += 1
        return self

    def __mul__(self, other):
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
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Multiplication with type {type(other)} is undefined.")

        r = copy.deepcopy(self.diffpat)

        for i in range(len(r)):
            r[i] = r[i] * other
        return DiffractionPattern(r, do_deep_copy=False)

    def __truediv__(self, other):
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
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        r = copy.deepcopy(self.diffpat)

        for i in range(len(r)):
            r[i] = r[i] / other
        return DiffractionPattern(r, do_deep_copy=False)

    def __floordiv__(self, other):
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
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        r = copy.deepcopy(self.diffpat)

        for i in range(len(r)):
            r[i] = r[i] // other
        return DiffractionPattern(r, do_deep_copy=False)

    def __add__(self, other):
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
        if isinstance(other, float) or isinstance(other, int):
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] + other
            return DiffractionPattern(r, do_deep_copy=False)
        elif isinstance(other, list):
            if len(other) != len(self.diffpat):
                raise ValueError(
                    f"You gave a list of {len(other)} to subtract from {len(self)} intensities. They need to be the same.")
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] + other[i]
            return DiffractionPattern(r, do_deep_copy=False)
        else:
            return self._add_and(other, operator.add)

    def __sub__(self, other):
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
        if isinstance(other, float) or isinstance(other, int):
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] - other
            return DiffractionPattern(r, do_deep_copy=False)
        elif isinstance(other, list):
            if len(other) != len(self.diffpat):
                raise ValueError(
                    f"You gave a list of {len(other)} to subtract from {len(self.diffpat)} intensities. They need to be the same.")
            r = copy.deepcopy(self.diffpat)
            for i in range(len(r)):
                r[i] = r[i] - other[i]
            return DiffractionPattern(r, do_deep_copy=False)
        else:
            return self._add_and(other, operator.sub)

    def __and__(self, other):
        """
        Overriding & operator to mean averagein diffraction pattern intensities together if the same angle

        Parameters
        ----------
        other : DiffractionPattern.

        Returns
        -------
        DiffractionPattern
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._add_and(other, operator.and_)

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
        if isinstance(other, float) or isinstance(other, int):
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
        if isinstance(other, float) or isinstance(other, int):
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

    def __iand__(self, other):
        """
        Overriding &= operator to mean averaging diffraction pattern intensities together if the same angle

        Parameters
        ----------
        other : DiffractionPattern.

        Returns
        -------
        DiffractionPattern
            with the intensities added together, and errors done nicely.

        Raises
        ------
        ValueError if wrong type or angles not equal used.

        """
        return self._add_and(other, operator.iand)

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
        if not (isinstance(other, float) or isinstance(other, int)):
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
        if not (isinstance(other, float) or isinstance(other, int)):
            raise TypeError(f"Division with type {type(other)} is undefined.")

        for i in range(len(self.diffpat)):
            self.diffpat[i] //= other
        return self

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class InterpolatedDiffractionPattern(DiffractionPattern):

    def __init__(self, data, interp_stepsize, aveStepSize=None, do_deep_copy=True):
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

        super().__init__(data, aveStepSize, do_deep_copy)
        # self.orig_diffpat = copy.deepcopy(self.diffpat)
        self.diffpat = InterpolatedDiffractionPattern.spline(interp_stepsize, self)

        # update ave step size after splining
        n = 0
        tot = 0
        for i in range(1, len(self.diffpat)):
            a = self.diffpat[i - 1].angle
            b = self.diffpat[i].angle
            n += 1
            tot += b - a
        self.aveStepSize = tot / n

    # def getOrigData(self):
    #     return copy.deepcopy(self.orig_diffpat)

    def spline(step_spline, dp):
        # print(f"Now splining {dp.filename} in steps of {step_spline}")
        IDP = InterpolatedDiffractionPattern
        # these will hold the full splined data once I've finished
        angle_spline = []
        intensity_spline = []
        error_spline = []

        # These are the original data that I am going to spline
        dp_angle = dp.getAngles()
        dp_intensity = dp.getIntensities()
        dp_error = dp.getErrors()

        # arbitrary choice of cut-off size -  If i encounter a stepsize
        #   greater then 5*the average, then I'll assume that is the next module
        step_threshold = 5 * dp.aveStepSize

        # Now I start the spline process
        start_index = 0
        stop_index = 0
        for i in range(1, len(dp_angle)):
            step = dp_angle[i] - dp_angle[i - 1]

            if step < step_threshold and i != len(dp_angle) - 1:
                stop_index += 1
            else:
                stop_index += 1  # to make stop_index exclusive

                # this is the range of data I want to interpolate
                #  I want to trim off a few datapoints either side of the module edge
                #  so I don't have to deal with their noisy edges.
                DROP_POINTS = 5
                angle_list = dp_angle[start_index + DROP_POINTS:stop_index - DROP_POINTS]
                intensity_list = dp_intensity[start_index + DROP_POINTS:stop_index - DROP_POINTS]
                error_list = dp_error[start_index + DROP_POINTS:stop_index - DROP_POINTS]

                # this is the interpolated data
                angle_interp = IDP.generateInterpList(angle_list[0], angle_list[-1], step_spline)
                intensity_interp = IDP.cubic_interp1d(angle_interp, angle_list, intensity_list)
                error_interp = IDP.cubic_interp1d(angle_interp, angle_list, error_list)

                # put the just-interpolated-data into the whole list of data to be used to make the new DP
                angle_spline += angle_interp
                intensity_spline += intensity_interp
                error_spline += error_interp

                start_index = stop_index

        # assemble it all into a XDP list
        r = []
        for i in range(len(angle_spline)):
            r.append(XRayDataPoint(angle_spline[i], intensity_spline[i], error_spline[i]))

        return r

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
        r = []

        angle = dp.getAngles()
        intensity = dp.getIntensities()
        error = dp.getErrors()

        # check limits on interpolation
        if interp[0] < angle[0]:  # ie I'm interpolating before the data starts
            raise ValueError(f"Trying to interpolate out of range: {interp[0]} < {angle[0]}")
        elif interp[-1] > angle[-1]:  # ie I'm interpolating after the data ends
            raise ValueError(f"Trying to interpolate out of range: {interp[-1]} > {angle[-1]}")
        else:
            pass  # It's all good, mate!

        # interpolate intensity and error
        intensity_interp = InterpolatedDiffractionPattern.cubic_interp1d(interp, angle, intensity)
        error_interp = InterpolatedDiffractionPattern.cubic_interp1d(interp, angle, error)

        # construct the diffpat
        for i in range(len(interp)):
            r.append(XRayDataPoint(interp[i], intensity_interp[i], error_interp[i]))

        return r

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

    def cubic_interp1d(x0, x, y, do_checks=True):
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

        def diff(lst):
            """
            numpy.diff with default settings
            """
            size = len(lst) - 1
            r = [0] * size
            for i in range(size):
                r[i] = lst[i + 1] - lst[i]
            return r

        def list_searchsorted(listToInsert, insertInto):
            """
            numpy.searchsorted with default settings
            """

            def float_searchsorted(floatToInsert, insertInto):
                """
                Helper function
                """
                for i in range(len(insertInto)):
                    if floatToInsert <= insertInto[i]:
                        return i
                return len(insertInto)

            return [float_searchsorted(i, insertInto) for i in listToInsert]

        def clip(lst, min_val, max_val, inPlace=False):
            """
            numpy.clip
            """
            if not inPlace:
                lst = lst[:]
            for i in range(len(lst)):
                if lst[i] < min_val:
                    lst[i] = min_val
                elif lst[i] > max_val:
                    lst[i] = max_val
            return lst

        def subtract(a, b):
            """
            returns a - b
            """
            return a - b

        if type(x0) is float:
            x0 = [x0]

        size = len(x)

        xdiff = diff(x)
        ydiff = diff(y)

        # allocate buffer matrices
        Li = [0] * size
        Li_1 = [0] * (size - 1)
        z = [0] * (size)

        # fill diagonals Li and Li-1 and solve [L][y] = [B]
        Li[0] = sqrt(2 * xdiff[0])
        Li_1[0] = 0.0
        B0 = 0.0  # natural boundary
        z[0] = B0 / Li[0]

        for i in range(1, size - 1, 1):
            Li_1[i] = xdiff[i - 1] / Li[i - 1]
            Li[i] = sqrt(2 * (xdiff[i - 1] + xdiff[i]) - Li_1[i - 1] * Li_1[i - 1])
            Bi = 6 * (ydiff[i] / xdiff[i] - ydiff[i - 1] / xdiff[i - 1])
            z[i] = (Bi - Li_1[i - 1] * z[i - 1]) / Li[i]

        i = size - 1
        Li_1[i - 1] = xdiff[-1] / Li[i - 1]
        Li[i] = sqrt(2 * xdiff[-1] - Li_1[i - 1] * Li_1[i - 1])
        Bi = 0.0  # natural boundary
        z[i] = (Bi - Li_1[i - 1] * z[i - 1]) / Li[i]

        # solve [L.T][x] = [y]
        i = size - 1
        z[i] = z[i] / Li[i]
        for i in range(size - 2, -1, -1):
            z[i] = (z[i] - Li_1[i - 1] * z[i + 1]) / Li[i]

        # find index
        index = list_searchsorted(x0, x)
        index = clip(index, 1, size - 1)

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
