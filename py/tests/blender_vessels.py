# -*- coding: utf-8 -*-
''' inspired by
https://github.com/njanakiev/blender-scripting.git
'''
'''
call with
blender -b -P blender_vessels.py
not -b for background
    -P for python
'''
import bpy
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

#find blender installation path
print("bpy is at: %s" % os.path.dirname(bpy.data.filepath))
def createSphere(origin=(0, 0, 0)):
    # Create icosphere
    bpy.ops.mesh.primitive_ico_sphere_add(location=origin)
    obj = bpy.context.object
    return obj
def createVessel(a=np.array((0,0,0)), b=np.array((1,0,0)), r=10):
  #calculate vessel length
  direction_vector = b-a;
  vessel_length = np.linalg.norm(b-a);
  print("vessel length: %f" % vessel_length)
  #calcualte center of mass
  print(a)
  print(b)
  print(b-a)
  midpoint = a + 0.5*(b-a);
  print("vessel center of mass: %s" % midpoint)
  print("euler angels:")
  beta =np.arccos(direction_vector[2])
  sin_beta = np.sqrt(1-direction_vector[2]*direction_vector[2])
  alpha=np.arctan2(direction_vector[0],direction_vector[2])
  padj=np.sqrt(direction_vector[0]*direction_vector[0]+direction_vector[2]*direction_vector[2])
  gamma=np.arctan2(padj,direction_vector[1])
  
  angel_with_xz_plane = np.arccos(np.dot([0,1,0],(b-a))/vessel_length)
  angel_with_xy_plane = pi/2-np.arccos(np.dot([0,0,1],(b-a))/vessel_length)
  #angel_with_z = 0
  print("yz: %f, xz: %f, xy: %f" % (alpha, beta,gamma))
  #bpy.ops.mesh.primitive_cylinder_add(radius=r, depth=2, end_fill_type='NOTHING', view_align=False, enter_editmode=False, location=(midpoint[0],midpoint[1],midpoint[2]), rotation=(alpha, beta, gamma))
  
  bpy.ops.mesh.primitive_cylinder_add(radius=r, depth=1, end_fill_type='NOTHING', view_align=False, enter_editmode=False, location=(midpoint[0],midpoint[1],midpoint[2]))
  #bpy.ops.mesh.primitive_cylinder_add(radius=r, depth=1, end_fill_type='NOTHING', view_align=False, enter_editmode=False, location=(0,0,0))
  bpy.ops.transform.resize(value=(1.0, 1.0, vessel_length), constraint_axis=(False, False, True), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
  bpy.context.object.rotation_mode = 'XYZ'
  bpy.context.object.rotation_euler[0] = alpha
  bpy.context.object.rotation_euler[1] = beta
  bpy.context.object.rotation_euler[2] = -gamma

  #bpy.ops.transform.rotate(value=pi/4, axis=(1, 0, 0), constraint_axis=(True, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
  #bpy.ops.transform.resize(value=(1.0, 1.0, vessel_length), constraint_axis=(False, False, True), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
  #bpy.ops.transform.rotate(value=angel_with_z, axis=(0, 0, 1), constraint_axis=(False, False, True), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
  #bpy.ops.transform.rotate(value=angel_with_x, axis=(1, 0, 0), constraint_axis=(True, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
  #bpy.ops.transform.rotate(value=angel_with_y, axis=(0, 1, 0), constraint_axis=(False, True, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
#  bpy.ops.transform.rotate(value=angel_with_z, axis=(0, 0, 1), constraint_axis=(True, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1, release_confirm=True, use_accurate=False)
  bpy.ops.object.convert(target='MESH')
  bpy.ops.object.modifier_add(type='SUBSURF')
  bpy.context.object.modifiers["Subsurf"].levels = 4
  bpy.context.object.modifiers["Subsurf"].render_levels = 4
  bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Subsurf")
  bpy.ops.object.editmode_toggle()
  bpy.ops.mesh.flip_normals()
  bpy.ops.object.editmode_toggle()
  
  
  #create displacement texture by clouds
  displacetex = bpy.data.textures.new('Texture', 'CLOUDS')
  displacetex.noise_scale = 2.1
  
  # Create and apply displace modifier
  displace = bpy.context.object.modifiers.new('Displace', 'DISPLACE')
  displace.texture = displacetex
  #displace.texture_coords = 'OBJECT'
  #displace.texture_coords_object = empty
  displace.mid_level = 0.5
  displace.strength = 3.6

  obj = bpy.context.object
  return obj

def createVessel_Path():
  bpy.ops.curve.primitive_nurbs_path_add(radius=1, view_align=False, enter_editmode=False, location=(0, 0, 0))

if __name__ == '__main__':
    # Remove all elements
    utils.removeAll()

    # Create camera
    bpy.ops.object.add(type='CAMERA', location=(0, 0, 20))
    cam = bpy.context.object
    cam.rotation_euler = Euler((0, 0, pi/2), 'XYZ')
    # Make this the current camera
    bpy.context.scene.camera = cam

    # Create lamps
    utils.rainbowLights()

    # Create object and its material
    #sphere = createSphere()
    #utils.setSmooth(sphere, 3)
    vessel_a = createVessel(a=np.array((0,0,0)), b=np.array((10,10,0)), r=2.5)
    vessel_b = createVessel(a=np.array((0,0,0)), b=np.array((10,-10,0)), r=2.5)
    vessel_c = createVessel(a=np.array((0,0,0)), b=np.array((-10,0,0)), r=3.15)
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
