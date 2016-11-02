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
#ifndef MY_ASSERT_H
#define MY_ASSERT_H

// my own Assert
#ifdef assert
#  undef assert
#endif
void _Assert( const char *msg, const char* file, int line );
#ifdef DEBUG
#define assert(exp) (void)( (exp) || (_Assert(#exp, __FILE__, __LINE__), 0) )
#define myAssert(exp) (void)( (exp) || (_Assert(#exp, __FILE__, __LINE__), 0) )
#define AssertMsg(str) (_Assert(str,__FILE__,__LINE__))
#else
#define assert(exp) ((void)0)
#define myAssert(exp) ((void)0)
#define AssertMsg(str) ((void)0)
#endif

#define AssertAlways(exp) (void)( (exp) || (_Assert(#exp, __FILE__, __LINE__), 0) )
#define AssertAlwaysMsg(msg) (_Assert(msg,__FILE__,__LINE__))


#endif