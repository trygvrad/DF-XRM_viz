import three_d_draw_object_classes as object_classes
import three_d_draw
import numpy as np

def add_crystal_structure( cif_file, scale = 4.0, rotation_function = None, legend_pos_shift = [0,0,0], cage_list = ['Ti', 'Nb'], axes_shift=None, make_bonds = ['X','X'], linewidth=1,
                            max_bond_length = 2.5, min_bond_length = 0,
                            oxygen_cage_radius = 2.5,
                            bounding_box_facecolor = [0,0,0,0],
                            cage_line_color = [0,0,0,1],
                            linecolor = [1,1,1,1],
                            show_text = True,
                            atom_legend_step = np.array([1,0,0])):
    '''
    adds a crystal structure
    input:
        cif_file: file path for the cif file
        scale: scale of the crystal structure
        rotation_function: a function rotation_function() that accepts a vector [x,y,z] and returns the rotated vector [xr,yr,zr]
                            This function describes the relation between the sample geometry and the crrystal geometry in the cif file
        legend_pos_shift: optional vector [x,y,z] that shifts the location of the legend
        cage_list: list of ions to make oxygen cages around
        axes_shift: default None, position of the axes of form [0,2,0] etc.
        make_bonds: atoms to make bonds between, i.e. ['Ti','O']
        show_text: if True, show text
        atom_legend_step: step for each atom in the legend, should be a vector pointing up in the final figure
    returns:
        list of three_d_draw_object_classes objects
    '''
    import Dans_Diffraction
    xtl = Dans_Diffraction.Crystal(cif_file)
    xtl.generate_lattice()
    # get atom cords and types
    loc_in_lattice_coord, atom_types, label, occ, uiso, mxmymz = xtl.Structure.generate_lattice(1, 1, 1) # get 3x3x3 unit cells of atoms
    # get the atoms we want to plot
    atoms_inside = np.sum((loc_in_lattice_coord>=-0.1)*(loc_in_lattice_coord<=1.1), axis=1)==3 # atoms inside single unit cell
    loc_in_lattice_coord = loc_in_lattice_coord[atoms_inside]
    atom_types = atom_types[atoms_inside]
    # center on zero
    loc_in_lattice_coord -= 0.5
    # calculate positions [in AA] in carthesian coordinates
    loc_in_AA = xtl.Cell.calculateR(loc_in_lattice_coord)
    # make atom objects
    objects = []
    for a_loc, a_typ in zip(loc_in_AA, atom_types):
        atom = object_classes.Atom(a_loc, a_typ)
        atom*=scale
        objects.append(atom)
    ############## make unit cell box #############
    node_locations=np.array([
        [-0.5,-0.5,-0.5],
        [ 0.5,-0.5,-0.5],
        [ 0.5, 0.5,-0.5],
        [-0.5, 0.5,-0.5],
        [-0.5,-0.5, 0.5],
        [ 0.5,-0.5, 0.5],
        [ 0.5, 0.5, 0.5],
        [-0.5, 0.5, 0.5],
    ])
    node_locations = xtl.Cell.calculateR(node_locations)
    node_locations *= scale
    # add box nodes
    box_node_collection = object_classes.NodeCollection(node_locations)
    # add box faces
    '''if len(facecolor) < 4 or facecolor[3] > 0:
        for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
        objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = facecolor, abscolor=abscolor))'''
    # add box lines
    for seg in [ [0,1], [1,2], [2,3], [3,0], [4,5], [5,6],[6,7], [7,4], [0,4], [1,5], [2,6], [3,7] ]:
        objects.append(object_classes.BoxLine(box_node_collection, seg, meshcolor = linecolor, linewidth = linewidth))

    if len(bounding_box_facecolor) < 4 or bounding_box_facecolor[3] > 0:
        for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
            objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = bounding_box_facecolor, abscolor=False))
    ############## add bonds ##############
    for i, typ in enumerate(atom_types):
        if typ == make_bonds[0]:
            for j, typ2 in enumerate(atom_types):
                if i == j:
                    continue
                #j = j+i+1
                if typ2 == make_bonds[1]:
                    loc0 = loc_in_AA[i]
                    loc1 = loc_in_AA[j]
                    #print(loc0,loc1,np.sqrt(np.sum((loc0-loc1)**2)) )
                    if np.sqrt(np.sum((loc0-loc1)**2)) < max_bond_length and np.sqrt(np.sum((loc0-loc1)**2)) > min_bond_length:
                        node_collection = object_classes.NodeCollection([loc0, loc1])
                        node_collection*=scale
                        objects.append(object_classes.BoxLine(node_collection, [0,1], meshcolor = [0.5,0.5,0.5,1], linewidth = linewidth))

    ############## add oxygen cages #############

    for center_loc, typ in zip(loc_in_AA, atom_types):
        if typ in cage_list:
            cage_in_AA = []
            for edge_loc, typ in zip(loc_in_AA, atom_types):
                if typ == 'O':
                    if np.sqrt(np.sum((center_loc-edge_loc)**2)) < oxygen_cage_radius:
                        cage_in_AA.append(edge_loc)
            node_collection = object_classes.NodeCollection(cage_in_AA)
            #print(np.array(cage_in_AA))
            #print(len(np.array(cage_atom_locs)))
            faces = []
            for i, loc_0 in enumerate(cage_in_AA):
                for j, loc_1 in enumerate(cage_in_AA[i+1:]):
                    vec = loc_1-loc_0
                    ortho = center_loc -0.5*(loc_1+loc_0)
                    if np.sqrt(np.sum(ortho**2))>1 and np.sqrt(np.sum(vec**2))<4:
                        #print([i,i+1+j],np.sqrt(np.sum(ortho**2)),np.sqrt(np.sum(vec**2)))
                        objects.append(object_classes.BoxLine(node_collection, [i,i+1+j], meshcolor = cage_line_color, linewidth = linewidth))
                        for k, loc_2 in enumerate(cage_in_AA[i+1+j+1:]):
                            vec0 = loc_2-loc_0
                            ortho0 = center_loc -0.5*(loc_2+loc_0)
                            vec1 = loc_2-loc_1
                            ortho1 = center_loc -0.5*(loc_2+loc_1)
                            if np.sqrt(np.sum(ortho0**2))>1 and np.sqrt(np.sum(vec0**2))<4:
                                if np.sqrt(np.sum(ortho1**2))>1 and np.sqrt(np.sum(vec1**2))<4:
                                    faces.append([i,i+1+j,i+1+j+1+k])
            node_collection*=scale
            if len(faces)>0:
                faces=np.array(faces)
                mesh = object_classes.Mesh(cage_in_AA, faces, meshcolor=[0,0,0,0], facecolor=[0.3,0.5,1,0.1])
                objects.append(mesh)
                mesh*=scale
                #print(np.array(faces))
    ############### make axes #############
    origo = np.array([0,0,0])
    axes_objects=[]
    # 100
    end = xtl.Cell.calculateR([1,0,0])[0]
    end = end/np.sqrt(np.sum(end**2))*2.0*scale
    mesh = three_d_draw.make_arrow_mesh(origo, origo+end, head_fraction = 0.45, rod_radius = 0.1, hat_radius = 0.25, nodes_per_circle = 12, color=[1,0,0,1], angle_offset=0)
    axes_objects.append(mesh)
    if show_text:
        text = object_classes.Text(origo+end*1.2, 'a', color=[0,0,0,1], scale=0.9)
        axes_objects.append(text)
    # 010
    end = xtl.Cell.calculateR([0,1,0])[0]
    end = end/np.sqrt(np.sum(end**2))*2.0*scale
    mesh = three_d_draw.make_arrow_mesh(origo, origo+end, head_fraction = 0.45, rod_radius = 0.1, hat_radius = 0.25, nodes_per_circle = 12, color=[0,1,0,1], angle_offset=0)
    axes_objects.append(mesh)
    if show_text:
        text = object_classes.Text(origo+end*1.2, 'b', color=[0,0,0,1], scale=0.9)
        axes_objects.append(text)
    # 010
    end = xtl.Cell.calculateR([0,0,1])[0]
    end = end/np.sqrt(np.sum(end**2))*2.0*scale
    mesh = three_d_draw.make_arrow_mesh(origo, origo+end, head_fraction = 0.45, rod_radius = 0.1, hat_radius = 0.25, nodes_per_circle = 12, color=[0,0,1,1], angle_offset=0)
    axes_objects.append(mesh)
    text = object_classes.Text(origo+end*1.2, 'c', color=[0,0,0,1], scale=0.9)
    if show_text:
        axes_objects.append(text)
        origo = np.min(box_node_collection.node_locations, axis=0)*1.1
    for obj in axes_objects:
        objects.append(obj)
    ############## rotate #############
    if not type(rotation_function) == type(None):
        node_collections=[]
        for obj in objects:
            if hasattr(obj,'node_collection'):
                if not obj.node_collection in node_collections:
                    node_collections.append( obj.node_collection )
                    nodes = obj.node_collection.node_locations
                else:
                    continue
            elif hasattr(obj,'node_locations'):
                    nodes = obj.node_locations
            else:
                nodes = [obj.loc]
            for node in nodes:
                rotation_function(node)
    ############## make legend #############
    legend_pos = np.array([ np.max(box_node_collection.node_locations[:,0]),
                            np.min(box_node_collection.node_locations[:,1]),
                            np.min(box_node_collection.node_locations[:,2])])
    legend_pos = legend_pos+np.array(legend_pos_shift)*scale
    legend_atoms = []
    # sort accrding to position in name of cif_file
    atom_type_set = set(atom_types)
    atom_type_order = []
    for atom in atom_type_set:
        atom_type_order.append(cif_file.split('/')[-1].find(atom))
    sorted_types = sorted(zip(atom_type_order,atom_type_set))
    atom_type_order, atom_type_set = [ list(tuple) for tuple in zip(*sorted_types)]

    for i, a_typ in enumerate(atom_type_set):
        atom = object_classes.Atom(np.array([0.0,0.0,0.0]), a_typ)
        objects.append(atom)
        atom *= scale
        atom += legend_pos
        atom += -2*atom_legend_step*i*scale
        if show_text:
            text = object_classes.Text(np.copy(atom.loc)+np.array([0,0,0.8])*scale, a_typ, color=[0,0,0,1], scale=0.9)
            objects.append(text)
            #text+=np.array([-2*i,0,0])*scale
    ############## shift axes #############
    if type(axes_shift) == type(None):
        origo = np.max(box_node_collection.node_locations, axis=0)*1.1+np.array([0,0,2])*scale
        origo[0] = np.min(box_node_collection.node_locations, axis=0)[0]*1.1
        for obj in axes_objects:
            obj += origo+np.array([0,2,0])*scale
    else:
        for obj in axes_objects:
            obj += np.array(axes_shift)*scale
    return objects
