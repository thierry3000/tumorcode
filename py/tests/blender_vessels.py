# -*- coding: utf-8 -*-
''' inspired by
https://github.com/njanakiev/blender-scripting.git
https://www.dailymotion.com/video/x2pcw13
'''
'''
call with
blender -b -P blender_vessels.py
not -b for background
    -P for python
'''
import bpy
import bmesh
import numpy as np
from math import pi
from mathutils import Euler
tau = 2*pi

# Check if script is opened in Blender program
import os, sys
if(bpy.context.space_data == None):
    cwd = os.path.dirname(os.path.abspath(__file__))
else:
    cwd = os.path.dirname(bpy.context.space_data.text.filepath)
# Get folder of script and add current working directory to path
sys.path.append(cwd)
import utils

''' in tumor code each segment has a fix lenght.
for visual representation each segment will 
be splittet into 
num_seg
'''
num_seg=30
#find blender installation path
print("bpy is at: %s" % os.path.dirname(bpy.data.filepath))
def createSphere(origin=(0, 0, 0)):
    # Create icosphere
    bpy.ops.mesh.primitive_ico_sphere_add(location=origin)
    obj = bpy.context.object
    return obj


def createVessel_Path(a=np.array((0,0,0)), b=np.array((1,0,0)), r=2.5):
  #create simple circle template
  
  bpy.ops.curve.primitive_bezier_circle_add(radius=r,view_align=False, enter_editmode=False, location=(0, 0, 0))
  dummy_circ = bpy.context.object
  print("dummy_circ:")
  print(dummy_circ)
  #curve to draw the circle along
  cu = bpy.data.curves.new("MyCurveData", "CURVE")
  
  pseudoCylinder = bpy.data.objects.new("MyCurveObject", cu)
  polyline = cu.splines.new('POLY')

  polyline.points.add(num_seg-1)
  '''calcultion will be needed here
  from point a to point b
  '''
  diff_vec = b-a;
  for i in range(num_seg):
    this_point = a+(i/num_seg)*diff_vec;
    polyline.points[i].co = (this_point[0],this_point[1],this_point[2],1)

  cu.dimensions = '3D'
  cu.bevel_object = bpy.data.objects["BezierCircle"]
  
  #bpy.context.scene.objects.link(pseudoCylinder)
  #bpy.context.scene.objects.active= pseudoCylinder
  
  
  me = pseudoCylinder.to_mesh(bpy.context.scene,False, 'PREVIEW')
  # add an object
  o = bpy.data.objects.new("MBallMesh", me)
  bpy.context.scene.objects.link(o)
  o.matrix_world = pseudoCylinder.matrix_world

    # not keep original
  #bpy.context.scene.objects.unlink(o)
  
  
  bpy.context.scene.objects.active = o
  
  #decimate
#  decimate = bpy.context.object.modifiers.new('Decimate', 'DECIMATE')
#  print(decimate)
#  decimate.decimate_type='UNSUBDIV'
#  decimate.iterations=0
#  bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Decimate")
  
  #subsurface
#  bpy.ops.object.modifier_add(type='SUBSURF')
#  bpy.context.object.modifiers["Subsurf"].levels = 2
#  bpy.context.object.modifiers["Subsurf"].render_levels = 2
#  bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Subsurf")
  
  #flip normals
#  bpy.ops.object.editmode_toggle()
#  bpy.ops.mesh.flip_normals()
#  bpy.ops.object.editmode_toggle()
  
  
  #create displacement texture by clouds
  displacetex = bpy.data.textures.new('Texture', 'CLOUDS')
  displacetex.noise_scale = 1.1
  
  # Create and apply displace modifier
  displace = bpy.context.object.modifiers.new('Displace', 'DISPLACE')
  displace.texture = displacetex
  displace.texture_coords = 'OBJECT'
  #displace.texture_coords_object = empty
  displace.mid_level = 3.
  displace.strength = 6.
  bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Displace")
  
  obj = bpy.context.object
  return obj
  
if __name__ == '__main__':
    # Remove all elements
    utils.removeAll()

    bpy.ops.object.lamp_add(type='SUN', radius=1, view_align=False, location=(0, 0, 50), rotation=(0, 0, -1.5708), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))

    # Create camera
    bpy.ops.object.add(type='CAMERA', location=(35, 0, 80))
    cam = bpy.context.object
    cam.rotation_euler = Euler((0.3, 0, pi/2), 'XYZ')
    # Make this the current camera
    bpy.context.scene.camera = cam

    # Create lamps
    #utils.rainbowLights()

    # Create object and its material
    #sphere = createSphere()
    #utils.setSmooth(sphere, 3)
    vessel_a = createVessel_Path(a=np.array((0,0,0)), b=np.array((100,100,0)), r=2.5)
    vessel_b = createVessel_Path(a=np.array((0,0,0)), b=np.array((100,-100,0)), r=2.5)
    vessel_b = createVessel_Path(a=np.array((0,0,0)), b=np.array((-100,0,0)), r=3.15)
#    vessel_c = createVessel(a=np.array((0,0,0)), b=np.array((-10,0,0)), r=3.15)
#    vessel_d = createVessel(a=np.array((-20,0,0)), b=np.array((-10,0,0)), r=3.5)
    
    #utils.setSmooth(vessel,3)

    # Specify folder to save rendering
    render_folder = os.path.join(cwd, 'rendering')
    if(not os.path.exists(render_folder)):
        os.mkdir(render_folder)

    # Render image
    rnd = bpy.data.scenes['Scene'].render
    rnd.resolution_x = 500
    rnd.resolution_y = 500
    rnd.resolution_percentage = 100
    rnd.filepath = os.path.join(render_folder, 'simple_sphere.png')
    bpy.ops.render.render(write_still=True)
