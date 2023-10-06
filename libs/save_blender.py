
import bpy
import math
from mathutils import Vector, Matrix
import numpy as np
import inspect, os
import time





def create_quadrilateral(nodes, name):
    '''
    creates a quadilateral facet spanning the nodes in nodes
    input:
        nodes: numpy array of shape [4,3]
        name: string
    return:
        bpy.data.object, bpy.data.mesh
    '''
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
    '''
    line connecting the nodes in nodes
    input:
        nodes: numpy array of shape [2,3]
        name: string
        r: float, radius of the line, default 0.1
    return:
        bpy.data.object
    '''
    d10 = nodes[1]-nodes[0]
    dist = np.sqrt(np.sum(d10**2))

    mesh = cyl.copy()
    #bpy.ops.mesh.primitive_cylinder_add(
    #    radius = r,
    #    depth = dist,
    #    location = np.average(nodes,axis=0)
    #    )
    mesh.transform(Matrix.Diagonal([r, r, dist, 1]))
    mesh.name = name

    ob = bpy.data.objects.new(name, mesh)
    ob.location = np.average(nodes,axis=0)
    phi = math.atan2(d10[1], d10[0])
    theta = math.acos(d10[2]/dist)
    #myobject = bpy.context.object
    bpy.context.collection.objects.link(ob)
    ob.rotation_euler[1] = theta
    ob.rotation_euler[2] = phi
    return ob #myobject



def set_mat_keys(mat,color):
    '''
    sets material properties of material mat
    input:
        mat: blender material object
        color: color [r,g,b,a]
    return:
        None
    '''
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

def save_blend(objects, filepath):
    '''
    objects is a list of 3D objects stored as dictionaries

    '''
    bpy.ops.wm.open_mainfile(filepath="template.blend")
    main_collection = bpy.data.collections.new("DF_XRM_vis")
    bpy.context.scene.collection.children.link(main_collection)  # Link collection to scene

    # all rods are based on a single sylinder
    global cyl, collection
    bpy.ops.mesh.primitive_cylinder_add(
        vertices=32,
        radius=1,
        depth=1,
        enter_editmode=False,
    )
    cyl = bpy.context.object.data
    cyl.name = 'tmp'
    bpy.data.objects.remove(bpy.context.object, do_unlink=True)

    mat_dict_face = {}
    mat_dict_box = {}
    mat_dict_line = {}
    collection_dict = {}
    t0 = time.time()
    for i, object in enumerate(objects):
        if i%10 ==1:
            dt = time.time() - t0
            frac = i/len(objects)
            print(f'{i}, {frac*100:.3f} %, {dt:.1f} s, remaining {(1-frac)*dt/frac:.1f} s, {object["type"]}       end = \r')
        if object['type']=='Mesh':

            mesh = bpy.data.meshes.new(str(i))   # create a new mesh
            for f in mesh.polygons:
                    f.use_smooth = True
            ob = bpy.data.objects.new(str(i), mesh)      # create an object with that mesh
            ob.location = Vector((0,0,0)) #by.context.scene.cursor_location   # position object at 3d-cursor
            mat_str = str(object['facecolor'])
            if not mat_str in mat_dict_face.keys():
                mat = bpy.data.materials.new(str(i))
                ob.active_material = mat
                mat.diffuse_color = object['facecolor']
                mat.use_nodes=True
                nodes = mat.node_tree.nodes
                links = mat.node_tree.links
                mat_dict_face[mat_str] = mat
                set_mat_keys(mat,object['facecolor'])
            else:
                mat = mat_dict_face[mat_str]
            col_str = 'face'+mat_str+str(i//10000)
            if not col_str in collection_dict.keys():
                collection = bpy.data.collections.new("face_"+str(i))
                main_collection.children.link(collection)  # Link collection to scene
                collection_dict[col_str] = collection
            else:
                collection = collection_dict[col_str]

            collection.objects.link(ob)                # Link object to collection
            # st up material
            ob.active_material = mat
            nodes = np.stack((object['nodes'][:,0], -object['nodes'][:,2], object['nodes'][:,1] ), axis = -1)
            mesh.from_pydata(list(nodes*10**-3),[],list(object['faces']))   # edges or faces should be [], or you ask for problems
            mesh.update(calc_edges=True)    # Update mesh with new data

        elif object['type']=='BoxFacet':
            mat_str = str(object['facecolor'])
            if not mat_str in mat_dict_box.keys():
                mat = bpy.data.materials.new(str(i))
                mat.diffuse_color = object['facecolor']
                mat.use_nodes=True
                mat_dict_box[mat_str] = mat
                set_mat_keys(mat,object['facecolor'])
            else:
                mat = mat_dict_box[mat_str]
            col_str = 'box'+mat_str+str(i//10000)
            if not col_str in collection_dict.keys():
                collection = bpy.data.collections.new("box_"+str(i))
                main_collection.children.link(collection)  # Link collection to scene
                collection_dict[col_str] = collection
            else:
                collection = collection_dict[col_str]
            layer_collection = bpy.context.view_layer.layer_collection.children[main_collection.name].children[collection.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
            nodes = np.stack((object['span'][:,0], -object['span'][:,2], object['span'][:,1] ), axis = -1)
            ob, mesh = create_quadrilateral(nodes*10**-3, str(i))
            #bpy.context.collection.objects.unlink(ob)
            #collection.objects.link(ob)                # Link object to collection
            ob.active_material = mat

        elif object['type']=='BoxLine':
            #collection.objects.link(ob)                # Link object to collection
            #bpy.context.collection.objects.unlink(ob)
            mat_str = str(object['meshcolor'])
            if not mat_str in mat_dict_line.keys():
                mat = bpy.data.materials.new(str(i))
                mat.diffuse_color = object['meshcolor']
                mat.use_nodes=True
                mat_dict_line[mat_str] = mat
                set_mat_keys(mat,object['meshcolor'])
            else:
                mat = mat_dict_line[mat_str]
            col_str = 'line'+mat_str+str(i//10000)
            if not col_str in collection_dict.keys():
                collection = bpy.data.collections.new("lines_"+str(i))
                main_collection.children.link(collection)  # Link collection to scene
                collection_dict[col_str] = collection
            else:
                collection = collection_dict[col_str]

            layer_collection = bpy.context.view_layer.layer_collection.children[main_collection.name].children[collection.name]
            bpy.context.view_layer.active_layer_collection = layer_collection
            nodes = np.stack((object['span'][:,0], -object['span'][:,2], object['span'][:,1] ), axis = -1)
            ob = make_line(nodes*10**-3, str(i), r=object['linewidth']*10**-3)
            #collection.objects.link(ob)                # Link object to collection
            ob.active_material = mat

    else:
        ...#print(object['type'])
    bpy.data.meshes.remove(cyl)
    bpy.ops.wm.save_as_mainfile(filepath=filepath)
    print(f'saved {filepath}')
