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
#ifndef _ANY_TREE_H_
#define _ANY_TREE_H_
#include <boost/property_tree/ptree.hpp>
#include <boost/any.hpp>

namespace boost { namespace property_tree {


// from boost property_tree examples
template<class Ext, class Int = boost::any>
struct any_translator
{
    typedef Ext external_type;
    typedef Int internal_type;

    external_type get_value(const internal_type &value) const
    {
      return boost::any_cast<external_type>(value);
    }

    internal_type put_value(const external_type &value) const
    {
      return value;
    }
};

template<typename TheType>
struct translator_between<boost::any, TheType>
{
    typedef any_translator<TheType> type;
};


typedef basic_ptree<std::string, boost::any> atree;

ptree asptree(const atree &src);

}}

#endif