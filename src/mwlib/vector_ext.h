/**
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
*/
#ifndef VECTOR_EXT_H
#define VECTOR_EXT_H

#include <vector>
#include <cstddef>

template<class T>
inline T* get_ptr(std::vector<T> &v) { return v.size()>0 ? &v[0] : NULL; }

template<class T>
inline const T* get_ptr(const std::vector<T> &v) { return v.size()>0 ? &v[0] : NULL; }


#endif