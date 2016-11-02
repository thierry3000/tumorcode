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
#ifndef MY_LOG
#define MY_LOG

#include <iosfwd>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/categories.hpp>

#include <tbb/mutex.h>

#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <vector>

namespace my
{

namespace my_log_impl
{
  class MyLogSink
  {
    public:
      typedef std::string::value_type char_type;
      typedef boost::iostreams::sink_tag category;

      MyLogSink(std::ostream &os) : indent_next(true)
      {
        streambuff = os.rdbuf();
      }

      MyLogSink(const MyLogSink &other)
      {
        streambuff = other.streambuff;
        indent_dat = other.indent_dat;
        build_indent_string();
        indent_next = other.indent_next;
        // the mutex is not copied. it cannot. this is why there is this copy constructor.
      }

      void write_indent()
      {
        streambuff->sputn(indent_str.c_str(),indent_str.size());
        indent_next = false;
      }

      std::streamsize write(const char_type* s, std::streamsize n)
      {
        tbb::mutex::scoped_lock lock(mutex);
        std::streamsize num_written = 0;
        for (std::streamsize i = 0; i<n; ++i)
        {
          if (indent_next)
            write_indent();
          if (s[i] == '\n') { indent_next = true; }
          num_written += (s[i] == streambuff->sputc(s[i])) ? 1 : 0;
        }
        return num_written;
      }

      void push(const std::string &s)
      {
        tbb::mutex::scoped_lock  lock(mutex);
        indent_dat.push_back(s);
        build_indent_string();
      }

      void build_indent_string()
      {
        indent_str.clear();
        for (int i=0; i<indent_dat.size(); ++i)
        {
          indent_str.append(indent_dat[i]);
        }
      }
      
      void pop()
      {
        tbb::mutex::scoped_lock  lock(mutex);
        assert(!indent_dat.empty());
        indent_dat.pop_back();
        build_indent_string();
      }
      
    private:
      tbb::mutex mutex;
      std::vector<std::string> indent_dat;
      std::string indent_str;
      std::streambuf *streambuff;
      bool indent_next;
  };


  class Log : public boost::iostreams::stream<MyLogSink>
  {
      typedef boost::iostreams::stream<MyLogSink> Base;
    public:
      Log(std::ostream &os) : Base(os) { }
      Log() : Base(std::cout) {}
      void push(const std::string &s = "  ") { (*(Base*)this)->push(s); }
      void pop() { (*(Base*)this)->pop(); }
  };

  class LogScope : boost::noncopyable
  {
      Log &log;
    public:
      LogScope(Log &log_, const std::string &name) : log(log_) { log.push(name); }
      ~LogScope() { log.pop(); }
  };
}

using my_log_impl::Log;
using my_log_impl::LogScope;

}


#endif
