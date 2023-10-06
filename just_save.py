'''
Tile to test saving to a .blend file.
This functionality has to be in a separate process because streamlit is threaded
and bpy does not support threading.
But using multiprocessing supresses regular error messages.

This file exists to debug the .blend saving functionality while recieving
correct error messages.
'''


from libs import save_blender
import pickle

path = "DF_XRM_vis_to_blender.pickled"

objects = pickle.load( open( path, "rb" ) )
save_blender.save_blend(objects, 'export.blend')
