#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
# Copyright (c) 2008 Peter Shinners

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from operator import isNumberType
from operator import isSequenceType
from math import hypot,sqrt, sin, cos, atan2
getattr = getattr  # copy builtin to global for speed


__all__ = ["Vec2","Vec3"]


class Vec2(object):
    """Immutable floating point 2D vector class.
        The vector acts as an indexable
        and iterable pair of floating point numbers. Operators are
        overridden to perform elemntwise operations. The length related
        methods are safe against zero length vectors, and will not
        raise division exceptions. All vector arguments accept any number
        or object with a pair of numbers."""
        
    __slots__ = ("__vec2__",)

    def __init__(self, x, y=False):
        """Initialize a new vector from any number or pair of numbers."""
        if y is not False:
            if y is _InternalCall:
                self.__vec2__ = x
            else:
                self.__vec2__ = 0.0 + x, 0.0 + y
        elif isNumberType(x):
            x += 0.0
            self.__vec2__ = x, x
        else:
            x, y = getattr(x, "__vec2__", x)
            self.__vec2__ = 0.0 + x, 0.0 + y

    @property
    def x(self):
        return self.__vec2__[0]

    @property
    def y(self):
        return self.__vec2__[1]

    def __repr__(self):
        return "Vec2(%f, %f)" % self.__vec2__

    def __iter__(self):
        return self.__vec2__.__iter__()

    def __len__(self):
        """The number of dimensions for the vector, always two."""
        return 2

    def __getitem__(self, key):
        """Index the vector as a two item container"""
        return self.__vec2__[key]
    
    def length(self):
        """Length of the vector"""
        return hypot(*self.__vec2__)

    def length2(self):
        """Squared length of the vector"""
        vx, vy = self.__vec2__
        return vx * vx + vy * vy

    def setlength(self, value):
        """Assign the length of the vector. If the vector
            currently has a zero length vector, it is moved only
            in positive y."""
        vx, vy = self.__vec2__
        l = hypot(vx, vy)
        if l:
            l = value / l
            return Vec2((vx * l, vy * l), _InternalCall)
        return Vec2((0.0, l), _InternalCall)

    def normalize(self):
        """Set the length of the vector to exactly one, in place,
            unless the vector has zero length."""
        vx, vy = self.__vec2__
        l = hypot(vx, vy)
        if l:
            return Vec2((vx / l, vy / l), _InternalCall)
        return self
        
    def distance(self, other):
        """The distance between two vectors."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return hypot(vx - ox, vy - oy)
    
    def distance2(self, other):
        """The squared distance between two vectors."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return (vx - ox) ** 2.0 +  (vy - oy) ** 2.0

    def setdistance(self, other, distance):
        """Set the distance to another vector, in place. If the
            two vectors are in the same position, the vector is moved
            in positive y."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        dx = vx - ox
        dy = vy - oy
        l = hypot(dx, dy)
        if l:
            l = distance / l
            return Vec2((ox + vx * l, oy + vy * l), _InternalCall)
        return Vec2((ox, oy + distance), _InternalCall) 
 
    def perp(self):
        """Make the vector a perpendicular, and keep the same length.
            The perpendicular will be rotated 90 degrees clockwise."""
        vx, vy = self.__vec2__
        return Vec2((vy, -vx), _InternalCall)

    def rotate(self, degrees, pivot=None):
        radians = degrees * -0.017453292519943295  # pi/180
        c = cos(radians)
        s = sin(radians)
        vx, vy = self.__vec2__
        if pivot is None:
            return Vec2((vx*c - vy*s, vx*s + vy*c), _InternalCall)
        if isNumberType(pivot):
            ox = oy = pivot
        else:
            ox, oy = getattr(pivot, "__vec2__", pivot)
        vx -= ox
        vy -= oy
        x = vx*c - vy*s
        y = vx*s + vy*c
        return Vec2((x + ox, y + oy), _InternalCall)

    def angle(self, other=None):
        vx, vy = self.__vec2__
        if not vx or not vy:
            return 0.0
        if other is None:
            return atan2(vx, vy) * 57.295779513082323  # 180/pi
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        cross = vx * oy - vy * ox
        dot = vx * ox + vy * oy
        return atan2(cross, dot) * 57.295779513082323  # 180/pi

    def setangle(self, degrees, other=None):
        radians = degrees * 0.017453292519943295  # pi/180
        c = cos(radians)
        s = sin(radians)
        vx, vy = self.__vec2__
        if other is None:
            l = hypot(vx, vy)
            return Vec2((l*s, l*c), _InternalCall)
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        l = hypot(vx - ox, vy - oy)
        return Vec2((ox + l*s, oy + l*c), _InternalCall)

    def cross(self, other):
        """Cross multiplication with another vector."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return vx * oy - vy * ox

    def dot(self, other):
        """Dot product with another vector."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return vx * ox + vy * oy

    def atan2(self):
        """Pass the coordinates to the atan2 function"""
        return atan2(*self.__vec2__)

    def project(self, other):
        """Project onto another vector. The other vector does not
            need to be normalized.""" 
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
            length2 = other ** 2
            length2 += length2
        else:
            ox, oy = getattr(other, "__vec2__", other)
            length2 = ox ** 2 + oy ** 2
        if not length2:
            return self
        dot = vx * ox + vy * oy
        length = dot / length2
        return Vec2((ox * length, oy * length), _InternalCall)
            
    def reflect(self, other):
        """Reflect the direction of the vector based on the direction
            of another. The other vector does not need to be normalized."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
            length2 = other ** 2
            length2 += length2
        else:
            ox, oy = getattr(other, "__vec2__", other)
            length2 = ox ** 2 + oy ** 2
        if not length2:
            return self
        dot = vx * ox + vy * oy
        length = dot / length2 * 2
        return Vec2((ox * length - vx, oy * length - vy), _InternalCall)

    def almost(self, other, epsilon=1e-7):
        """Compare if two vectors are almost equal."""
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        dx = vx - ox
        dy = vy - oy
        nepsilon = -epsilon
        return (dx < epsilon and dy > nepsilon 
                    and dy < epsilon and dy > nepsilon)
        
    def floor(self):
        """Create a pair of integers of the floor from each value."""
        vx, vy = self.__vec2__
        _i = int
        return (_i(vx // 1.0), _i(vy // 1.0))  # NOTE, not a vec

    def ceil(self):
        """Create a pair of integers of the ceil from each value."""
        vx, vy = self.__vec2__
        _i = int
        return (_i(vx // 1.0) + 1, _i(vy // 1.0) + 1)  # NOTE, not a vec
        
    def __complex__(self):
        return complex(*self.__vec2__)
       
    def __neg__(self):
        vx, vy = self.__vec2__
        return Vec2((-vx, -vy), _InternalCall)

    def __pos__(self):
        return self

    def __abs__(self):
        vx, vy = self.__vec2__
        _a = abs
        return Vec2(_a(vx), _a(vy))

    def __invert__(self):
        vx, vy = self.__vec2__
        return Vec2((vy, -vx), _InternalCall)

    def __cmp__(self, other):
        if isNumberType(other):
            ovec = other, other
        else:
            ovec = getattr(other, "__vec2__", None)
            if ovec is None:
                ovec = tuple(other)                
        return cmp(self.__vec2__, ovec)

    def __nonzero__(self):
        vx, vy = self.__vec2__
        return bool(vx or vy)

    def __add__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((vx + ox, vy + oy), _InternalCall)

    __radd__ = __add__
        
    def __sub__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((vx - ox, vy - oy), _InternalCall)

    def __rsub__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((ox - vx, oy - vy), _InternalCall)
    
    def __mul__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((vx * ox, vy * oy), _InternalCall)
    
    __rmul__ = __mul__
        
    def __div__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((vx / ox, vy / oy), _InternalCall)

    def __rdiv__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((ox / vx, oy / vy), _InternalCall)
    
    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __floordiv__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((vx // ox, vy // oy), _InternalCall)

    def __rfloordiv__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((ox // vx, oy // vy), _InternalCall)
    
    __xor__ = dot
    __rxor__ = dot
    __mod__ = cross

    def __rmod__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return ox * vy - oy * vx

    def __pow__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            return Vec2((vx * other, vy * other), _InternalCall)
        else:
            ox, oy = getattr(other, "__vec2__", other)
            return Vec2((vx ** ox, vy ** oy), _InternalCall)

    def __rpow__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = other
        else:
            ox, oy = getattr(other, "__vec2__", other)
        return Vec2((ox ** vx, oy ** vy), _InternalCall)
    
    def __getstate__(self):
        return self.__vec2__
    
    def __setstate__(self, value):
        self.__vec2__ = value




class Vec3(object):
    """Immutable floating point 2D vector class.
        The vector acts as an indexable
        and iterable pair of floating point numbers. Operators are
        overridden to perform elemntwise operations. The length related
        methods are safe against zero length vectors, and will not
        raise division exceptions. All vector arguments accept any number
        or object with a pair of numbers."""
        
    __slots__ = ("__vec3__",)

    def __init__(self, x, y=None, z=None):
        """Initialize a new vector from any number or pair of numbers."""
        if y is _InternalCall:
            self.__vec3__ = x
        elif isinstance(x,Vec3):
            self.__vec3__ = getattr(x,"__vec3__",x)
        elif isSequenceType(x):
            self.__vec3__ = 0.0 + x[0], 0.0 + x[1], 0.0 + x[2]    
        elif isNumberType(x):
            x += 0.0
            if y is None:
                y = x
            else:
                y += 0.0
            if z is None:
                z = x
            else:
                z += 0.0
            self.__vec3__ = x,y,z
        else:
            raise RuntimeError('cannot init Vec3 with '+str(x))

    @property
    def x(self):
        return self.__vec3__[0]

    @property
    def y(self):
        return self.__vec3__[1]

    @property
    def z(self):
        return self.__vec3__[2]

    def __repr__(self):
        return "Vec3(%f, %f, %f)" % self.__vec3__

    def __iter__(self):
        return self.__vec3__.__iter__()

    def __len__(self):
        """The number of dimensions for the vector, always two."""
        return 3

    def __getitem__(self, key):
        """Index the vector as a two item container"""
        return self.__vec3__[key]
    
    def length(self):
        """Length of the vector"""
        return sqrt(self.length2())

    def length2(self):
        """Squared length of the vector"""
        vx, vy, vz = self.__vec3__
        return vx * vx + vy * vy + vz * vz

    def swapaxes(self, a, b):
        l = list(self.__vec3__)
        t = l[a]
        l[a] = l[b]
        l[b] = t
        return Vec3(tuple(l), _InternalCall)

    def replace(self, i, val):
        l = list(self.__vec3__)
        l[i] = val
        return Vec3(tuple(l),_InternalCall)

    def normalize(self):
        """Set the length of the vector to exactly one, in place,
            unless the vector has zero length."""
        vx, vy, vz = self.__vec3__
        l = self.length()
        if l:
            return Vec3((vx / l, vy / l, vz / l), _InternalCall)
        return self

    def cross(self, other):
        """Cross multiplication with another vector."""
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3(vy * oz - vz*oy ,vz*ox - vx*oz , vx * oy - vy * ox)

    def dot(self, other):
        """Dot product with another vector."""
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return vx * ox + vy * oy + vz * oz

    def almost(self, other, epsilon=1e-7):
        """Compare if two vectors are almost equal."""
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        dx = vx - ox
        dy = vy - oy
        dz = vz - oz
        nepsilon = -epsilon
        return (dx < epsilon and dy > nepsilon 
                    and dy < epsilon and dy > nepsilon)
        
    def floor(self):
        """Create a pair of integers of the floor from each value."""
        vx, vy, vz = self.__vec3__
        _i = int
        return (_i(vx // 1.0), _i(vy // 1.0), _i(vz // 1.0))  # NOTE, not a vec

    def ceil(self):
        """Create a pair of integers of the ceil from each value."""
        vx, vy, vz = self.__vec3__
        _i = int
        return (_i(vx // 1.0) + 1, _i(vy // 1.0) + 1, _i(vz // 1.0) + 1)  # NOTE, not a vec        
       
    def __neg__(self):
        vx, vy, vz = self.__vec3__
        return Vec3((-vx, -vy, -vz), _InternalCall)

    def __pos__(self):
        return self

    def __abs__(self):
        vx, vy, vz = self.__vec3__
        _a = abs
        return Vec3(_a(vx), _a(vy), _a(vz))

    def __cmp__(self, other):
        if isNumberType(other):
            ovec = other, other
        else:
            ovec = getattr(other, "__vec3__", None)
            if ovec is None:
                ovec = tuple(other)                
        return cmp(self.__vec3__, ovec)

    def __nonzero__(self):
        vx, vy, vz = self.__vec3__
        return bool(vx or vy or vz)

    def __add__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((vx + ox, vy + oy, vz + oz), _InternalCall)

    __radd__ = __add__
        
    def __sub__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((vx - ox, vy - oy, vz - oz), _InternalCall)

    def __rsub__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec2((ox - vx, oy - vy, oz - vz), _InternalCall)
    
    def __mul__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((vx * ox, vy * oy, vz * oz), _InternalCall)
    
    __rmul__ = __mul__
        
    def __div__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((vx / ox, vy / oy, vz / oz), _InternalCall)

    def __rdiv__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((ox / vx, oy / vy, oz / vz), _InternalCall)
    
    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __floordiv__(self, other):
        vx, vy = self.__vec2__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((vx // ox, vy // oy, vz // oz), _InternalCall)

    def __rfloordiv__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((ox // vx, oy // vy, oz // vz), _InternalCall)
    
    __xor__ = dot
    __rxor__ = dot
    __mod__ = cross

    def __rmod__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3(vy * oz - vz * oy ,vz * ox - vx * oz , vx * oy - vy * ox)

    def __pow__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            return Vec3((vx * other, vy * other), _InternalCall)
        else:
            ox, oy = getattr(other, "__vec2__", other)
            return Vec3((vx ** ox, vy ** oy, vz ** oz), _InternalCall)

    def __rpow__(self, other):
        vx, vy, vz = self.__vec3__
        if isNumberType(other):
            ox = oy = oz = other
        else:
            ox, oy, oz = getattr(other, "__vec3__", other)
        return Vec3((ox ** vx, oy ** vy, oz ** vz), _InternalCall)
    
    def __getstate__(self):
        return self.__vec3__
    
    def __setstate__(self, value):
        self.__vec3__ = value

    def max(self):
        return max(self.__vec3__)


class _InternalCall(object):
    """Passed to the Vec2 constructor for an optimized init case"""
        



