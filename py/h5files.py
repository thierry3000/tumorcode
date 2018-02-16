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
'''
    DEPRICATED SINCE h5py2.7
'''

## -*- coding: utf-8 -*-
#import h5py
#import atexit
#import os
#import sys
#
#if sys.flags.debug:
#    print("h5py.__file__: %s" % h5py.__file__)
#
#'''
#  Need this because - citing the h5py docs
#  http://www.h5py.org/docs/whatsnew/2.0.html#file-objects-must-be-manually-closed
#  "Please note that opening the same file multiple times (i.e. without closing it first) continues to result in undefined behavior.".
#  I need to access the same file multiple times though, but without knowing that it is in fact the same file!
#'''
#
#search_paths = map(lambda s: s.strip(), os.environ.get('H5FILES_SEARCH_DIRS', '') .split(':'))
#
#
#def FindFile_(filename, relatedH5Object):
#  dirs = []
#  if relatedH5Object is not None:
#    dirs.append(os.path.dirname(relatedH5Object.file.filename))
#  dirs += search_paths
#  for dir in dirs:
#    trial = os.path.join(dir, filename)
#    if os.path.isfile(trial):
#      return trial, None
#  return filename, Exception('%s attempt to find %s failed! Was looking in %s from %s! You can include search paths in the environment variable H5FILES_SEARCH_DIRS.' % (__name__, filename, dirs, os.getcwd()))
#
#
#class FileWrapper(h5py.File):
#  def __init__(self, fn, mode):
#    h5py.File.__init__(self, fn, mode)
#
#  def __enter__(self):
#    return self
#
#  def __exit__(self, type, value, traceback):
#    global close
#    close(self)
#
#
#class H5Files(object):
#  def __init__(self):
#    self.files = dict() # so we use a dict to keep track of the files
#
#  def makeKey_(self, fn):
#    _, ino, dev, _, _, _, _, _, _, _ = os.stat(fn)
#    return (ino, dev) # uniquely identifies a file on unix file systems
#    #http://effbot.org/zone/python-fileinfo.htm
#    #return os.path.realpath(os.path.abspath(os.path.normpath(fn)))
#
#  def open(self, fn, mode = 'r', relatedObjectForSearch = None, search = True):
#    if not os.path.isfile(fn) and (search or relatedObjectForSearch is not None):
#      fn, err = FindFile_(fn, relatedObjectForSearch)
#    else:
#      err = None
#    if os.path.isfile(fn): 
#      # does the file exist on disk, then generate a key and see if it is in the dict. 
#      #If not, then attempt to open it and put it in the dict after that.
#      key = self.makeKey_(fn)
#      try: 
#        f = self.files[key]
#        return f
#      except KeyError:
#        pass
#    # no luck, fallen through to here
#    try:
#      f = FileWrapper(fn, mode)
#    except Exception, e: # trouble opening the file
#      if err: 
#        print err.message
#      raise e
#    key = self.makeKey_(fn)    
#    self.files[key] = f
#    return f
#
#  def close(self, f):
#    key = self.makeKey_(f.filename)
#    del self.files[key]
#    f.close()
#
#  def closeall(self):
#    #for f in self.files.itervalues(): f.close()
#    self.files = dict()
#
#
#
#h5files_ = H5Files()
##obvious hack to close all files, is this still necessary with the new h5py?
#atexit.register(H5Files.closeall, h5files_)
#
#open = h5files_.open
#close = h5files_.close
#closeall = h5files_.closeall
#
#def openLink(g, name):
#  if g.id.links.get_info(name).type == h5py.h5l.TYPE_EXTERNAL:
#    filename, path = g.id.links.get_val(name)
#    f = h5files_.open(filename, g.file.mode, relatedObjectForSearch = g)
#    return f[path]
#  else:
#    return g[name]
#
#def isLink(g, name):
#  return True if  g.id.links.get_info(name).type in (h5py.h5l.TYPE_SOFT, h5py.h5l.TYPE_EXTERNAL) else False