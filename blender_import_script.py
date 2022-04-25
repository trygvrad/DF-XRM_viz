
import bpy
import pickle
import math
from mathutils import Vector
import numpy as np
import inspect, os

def dir_of_this():
    raw_path = inspect.getfile(inspect.currentframe())
    l = len(bpy.data.filepath)
    if raw_path[:l] == bpy.data.filepath:
        # Blender gives us a path consisting of the .blend file name plus the name of the text buffer inside the .blend.
        # Mildly awkward, but still usable.
        buf_name = raw_path[l+1:]
        #print("buffer is %s"%buf_name)
        raw_path = bpy.data.texts[buf_name].filepath
    return os.path.dirname(raw_path)


path = str(dir_of_this())+"/DF_XRM_vis_to_blender.pickled"

def create_quadrilateral(nodes, name):
    myvertex = []
    myfaces = []
    mypoint = [tuple(nodes[0])]
    myvertex.extend(mypoint)
    mypoint = [tuple(nodes[1])]
    myvertex.extend(mypoint)
    mypoint = [tuple(nodes[3])]
    myvertex.extend(mypoint)
    mypoint = [tuple(nodes[2])]
    myvertex.extend(mypoint)
    myface = [(0, 1, 3, 2)]
    myfaces.extend(myface)
    mymesh = bpy.data.meshes.new(name)
    myobject = bpy.data.objects.new(name, mymesh)
    collection.objects.link(myobject)
    # Generate mesh data
    mymesh.from_pydata(myvertex, [], myfaces)
    # Calculate the edges
    mymesh.update(calc_edges=True)
    # Set Location
    myobject.location.x = 0
    myobject.location.y = 0
    myobject.location.z = 0

    return myobject, mymesh

def make_line(nodes, name, r=0.1):
    d10 = nodes[1]-nodes[0]
    dist = np.sqrt(np.sum(d10**2))

    mesh = bpy.ops.mesh.primitive_cylinder_add(
        radius = r,
        depth = dist,
        location = np.average(nodes,axis=0)
        )


    phi = math.atan2(d10[1], d10[0])
    theta = math.acos(d10[2]/dist)
    myobject = bpy.context.object
    bpy.context.object.rotation_euler[1] = theta
    bpy.context.object.rotation_euler[2] = phi

    return myobject



def set_mat_keys(mat,color):
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    node_output  = nodes.new(type='ShaderNodeOutputMaterial')
    node_output.location = 400,0
    node_pbsdf = nodes.new(type='ShaderNodeBsdfPrincipled')
    node_pbsdf.location = 0,0
    node_pbsdf.inputs['Base Color'].default_value = color
    node_pbsdf.inputs['Alpha'].default_value = color[-1] # 1 is opaque, 0 is invisible
    node_pbsdf.inputs['Roughness'].default_value = 0.2
    node_pbsdf.inputs['Specular'].default_value = 0.5
    node_pbsdf.inputs['Transmission'].default_value = 0.5 # 1 is fully transparent

    link = links.new(node_pbsdf.outputs['BSDF'], node_output.inputs['Surface'])
    mat.node_tree.nodes["Principled BSDF"].inputs['Base Color'].default_value = color
    mat.node_tree.nodes["Principled BSDF"].inputs['Alpha'].default_value = color[-1]
    mat.node_tree.nodes["Principled BSDF"].inputs['Roughness'].default_value = 0.2
    mat.node_tree.nodes["Principled BSDF"].inputs['Specular'].default_value = 0.5
    mat.node_tree.nodes["Principled BSDF"].inputs['Transmission'].default_value = 0.5


objects = pickle.load( open( path, "rb" ) )

collection = bpy.data.collections.new("new_collection")
bpy.context.scene.collection.children.link(collection)                # Link object to scene

for i, object in enumerate(objects):
    if object['type']=='Mesh':

        mesh = bpy.data.meshes.new(str(i))   # create a new mesh
        for f in mesh.polygons:
                f.use_smooth = True
        ob = bpy.data.objects.new(str(i), mesh)      # create an object with that mesh
        ob.location = Vector((0,0,0)) #by.context.scene.cursor_location   # position object at 3d-cursor

        collection.objects.link(ob)                # Link object to collection
        # st up material
        mat = bpy.data.materials.new(str(i))
        ob.active_material = mat
        mat.diffuse_color = object['facecolor']
        mat.use_nodes=True
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        ob.active_material = mat
        set_mat_keys(mat,object['facecolor'])


        mesh.from_pydata(list(object['nodes']*10**-3),[],list(object['faces']))   # edges or faces should be [], or you ask for problems
        mesh.update(calc_edges=True)    # Update mesh with new data

    elif object['type']=='BoxFacet':
        ob, mesh = create_quadrilateral(object['span']*10**-3, str(i))
        mat = bpy.data.materials.new(str(i))
        mat.diffuse_color = object['facecolor']
        mat.use_nodes=True
        ob.active_material = mat
        set_mat_keys(mat,object['facecolor'])

    elif object['type']=='BoxLine':
        ob = make_line(object['span']*10**-3, str(i), r=object['linewidth']*10**-3)
        collection.objects.link(ob)                # Link object to collection
        bpy.context.collection.objects.unlink(ob)
        mat = bpy.data.materials.new(str(i))
        mat.diffuse_color = object['meshcolor']
        mat.use_nodes=True
        ob.active_material = mat
        set_mat_keys(mat,object['meshcolor'])

    else:
        print(object['type'])
