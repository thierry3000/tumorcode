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

#http://zulko.github.io/blog/2013/09/27/read-and-write-video-frames-in-python-using-ffmpeg/

''' usage:
  cat ../tum-q2d-mini-11layer-P10-typeA-sample00-q2d_11layer_off-out000*.png |
  python2 ~/tumorcode/py/create_movie.py -l typeA

1) create your images with usual povray
2) maybe merge adaption an no adaption
  for i in {1..5};
  do convert "input_image1${i}_.png" "input_image1${i}_.png" 
  +append "both_000${i}.png";done
  NOTE: THE BRACKETS
3) use this scritp see usage
'''

import subprocess as sp
import glob
import shlex
import datetime

import identifycluster
if(identifycluster.getname()=='snowden'):
  FFMPEG_BIN = "ffmpeg_latest"
else:
  FFMPEG_BIN = "ffmpeg" # on Linux ans Mac OS

if __name__ == '__main__':

  import optparse  #Note: Deprecated since version 2.7. Use argparse instead
  parser = optparse.OptionParser()
  parser.add_option("-l","--label", dest="label", help="attach label to output file", default='a_label')
  parser.add_option("-o","--overlay_merge", dest="overlay_merge",help="adds an overlay to distinguishe files", default = False, action="store_true")
#  parser.add_option("-T","--only_two_root", dest="two", help="flag to change the considered types", default=False, action="store_true")  
#  parser.add_option("-a","--with_all_types", dest="all_types", help="take all types",default=False, action="store_true")  
#  parser.add_option("-s","--singel_type", dest="single", help="", default=False, action="store_true")  
  options, args = parser.parse_args()  
  
  if 1:
    print("FFMPEG_BIN: %s" % FFMPEG_BIN)
    stamp = datetime.datetime.now().time()

    command = [ FFMPEG_BIN,
          '-y', # (optional) overwrite output file if it exists
          '-framerate', '1/0.3',
          '-i', '-', # The imput comes from a pipe        
          '-an', # Tells FFMPEG not to expect any audio
          '-c:v', 'libx264',
          '-r', '50',
          '-pix_fmt', 'yuv420p']    
    
    if( not options.label == 'a_label'):
      afilename = 'my_output_videofile_%s_%s.mp4' % (stamp, options.label)
    else:
      afilename = 'out_%s.mp4' % stamp
      
    if( options.overlay_merge):
      print("we are in overlay")
      command.append('-vf')
      command.append("drawtext=fontfile=/usr/share/fonts/libertine-ttf/LinBiolinum_aBL.ttf:text=\'Adaption\':fontcolor=white@0.9:fontsize=100:x=(w-text_w)/4:y=(h-text_h)/10")
      
    command.append(afilename)      
    print(command)
    
    sp.call(command)
'''works
cat both_000*.png |ffmpeg -framerate 1/0.3 -i - -an -c:v libx264 -r 50 -pix_fmt yuv420p -vf "drawtext=fontfile=/usr/share/fonts/libertine-ttf/LinBiolinum_aBL.ttf:text='Adaption':fontcolor=white@0.9:fontsize=100:x=(w-text_w)/4:y=(h-text_h)/10" -codec:a copy -y test.mp4
'''
