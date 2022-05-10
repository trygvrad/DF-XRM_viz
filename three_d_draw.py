import numpy as np
import copy
import three_d_draw_object_classes as object_classes
import time
import streamlit as st
import scipy.optimize


def make_3d_perspective(params, export_blender_filepath=None,
                # deals with rendering the figure
                crystal_rotation_function = None,
                scale = 0.1,
                tilt_to_vector = None,
                # deals with additional objects to render
                atom_list = None, crysta_structure_position = [-70, -10, 60],
                extend_beam = (1000,0), beam_opacity = 1,
                extend_imaging_system=(0,1000), imaging_system_lw = 1,
                draw_scattered_beam = True,
                lens_scale = 0,
                bounding_box_facecolor = [0,0,0,0.3],
                # stages and arrows
                draw_stage = False,
                draw_axes = True,
                draw_scattering_arrows = True,
                # misc
                arrow_nodes_per_circle=12,
                show_text = True,
                ):
    '''
    Makes a matplotlib figure illustrating the experiment
    input:
        params: dictionary specifying the geometry
        export_blender_filepath: string, optional filepath for blender export. Default None

        crystal_rotation_function: function tha takes an array of lenght three describing a vector in the crystal basis, and returns the vector in simulation coordiantes
        scale: float default = 0.1, multiplication factor for the simulated volume
        tilt_to_vector: vector of lenght 3 that defines 'up' or none

        atom_list: None or list of atoms and other objects to render crystal structure
        crysta_structure_position: position to render the crystal strcuture (with respect to the volume center of the simulated volume. This is in the coordinate system of the drawing: [horizontal, vertical, out-of-plane]
        extend_beam: tuple or None, default (1000,0). The beam is extended in (before, after) directions along the [001] axis equal to this value
                    If None, the beam is not shown
        beam_opacity: opacity of the beam
        extend_imaging_system: tuple (forward,back), default = (0,1000). The integrated field projection is extended in both directions along the [001] axis equal to this value.
                        It is invisible ulness imaging_system_lw is also set.
                        also determines the posistion of lenses, and detector_map
        imaging_system_lw: float, default is 0, linewith of the box of extended intensity
        draw_scattered_beam: bool default True, if true draws the scattered beam
        lens_scale: float, default=0, scale to render the lens in. If 0, do not show the lens
        bounding_box_facecolor: color of the sample bounding box, [r,g,b,a]

        draw_stage: bool default False. If True, draw the stage
        draw_axes: bool default True. If False, do not draw the arrows indicating the sample geometry
        draw_scattering_arrows: bool default True. If False, do not draw the arrows indicating the scattering vector

        arrow_nodes_per_circle: int default 12. Nodes per circel when drawing arrwos
        show_text: if True, show text
    returns:
        drawing, fig, ax
    '''
    outstr = []
    ax = {} # <- this list will be populated by plotly objects
    ax['data'] = []
    ax['annotations'] = []

    if not type(tilt_to_vector) == type(None):  # tilt fig rotation
        phi, chi, omega, beam_stage_angle = get_phi_chi_omega(params['k0_vector'], params['kh_vector'], tilt_to_vector )
    else:
        chi = 0
        phi = 0
        omega = 0
        beam_stage_angle = 0
    drawing = Drawing()
    # add crystal structure if cif specified
    if not type(atom_list) == type(None):
        crystal_structure = drawing.add_atom_list(atom_list)
    # add sample as box
    box_shape = params['shape']*scale
    sample_node_collection, _, __ = drawing.add_box(box_shape, [0,0,0], meshcolor=[0,0,0,1], linewidth = 2, facecolor=bounding_box_facecolor, abscolor = False)
    # text box size
    if show_text:
        drawing.add_text([0,1+0.5*box_shape[1],1+0.5*box_shape[2]], f'{box_shape[0]/scale:.2f} µm')
        drawing.add_text([1+0.5*box_shape[0],0,1+0.5*box_shape[2]], f'{box_shape[1]/scale:.2f} µm')
        drawing.add_text([3+0.5*box_shape[0],-0.5*box_shape[1],-8], f'{box_shape[2]/scale:.2f} µm')

    if draw_scattering_arrows:
        q_arrows_scale = 1.0
        Q_arrows = []
        Q_norm = params['Q_vector']/(np.sqrt(np.sum(params['k0_vector']**2))+np.sqrt(np.sum(params['kh_vector']**2)))
        ar = drawing.add_arrow([0,0,0], Q_norm*9*q_arrows_scale*100,color=[0,0,0,1],nodes_per_circle=arrow_nodes_per_circle,head_fraction = 0.45/1, rod_radius = 0.1, hat_radius = 0.25)
        Q_arrows.append(ar)
        if show_text:
            ar = drawing.add_text(Q_norm*9*q_arrows_scale*100, 'Q')
            Q_arrows.append(ar)
        k0_norm = params['k0_vector']/np.sqrt(np.sum(params['k0_vector']**2))
        ar = drawing.add_arrow(-k0_norm*4.5*q_arrows_scale*100, [0,0,0],color=[0.5,0,0,1],nodes_per_circle=arrow_nodes_per_circle,head_fraction = 0.45/3, rod_radius = 0.1/3, hat_radius = 0.25/3)
        k0_butt_node = ar.node_locations[0] # use this to rotate stage
        Q_arrows.append(ar)
        if show_text:
            ar = drawing.add_text((-k0_norm*4.5+np.array([1,0,0]))*q_arrows_scale*100, 'k0')
            Q_arrows.append(ar)
        kh_norm = params['kh_vector']/np.sqrt(np.sum(params['kh_vector']**2))
        ar = drawing.add_arrow([0,0,0], kh_norm*4.5*q_arrows_scale*100,color=[0.5,0,0,1],nodes_per_circle=arrow_nodes_per_circle,head_fraction = 0.45/3, rod_radius = 0.1/3, hat_radius = 0.25/3)
        Q_arrows.append(ar)
        if show_text:
            ar = drawing.add_text(kh_norm*4.5*q_arrows_scale*100, 'kh')
            Q_arrows.append(ar)

    # add id06 axes
    if draw_axes:
        xyz_arrows_scale = 1.0
        xyz_arrows = []
        ar = drawing.add_arrow([0,0,0], [0,0,3*xyz_arrows_scale*100],color=[1,0.1,0,1],nodes_per_circle=arrow_nodes_per_circle)
        xyz_arrows.append(ar)
        ar = drawing.add_arrow([0,0,0], [0,-3*xyz_arrows_scale*100,0],color=[0,0.9,0.3,1],nodes_per_circle=arrow_nodes_per_circle)
        xyz_arrows.append(ar)
        ar = drawing.add_arrow([0,0,0], [3*xyz_arrows_scale*100,0,0],color=[0,0.1,1,1],nodes_per_circle=arrow_nodes_per_circle, angle_offset=0)
        xyz_arrows.append(ar)
        if show_text:
            ar = drawing.add_text([3.2*xyz_arrows_scale*100,0,0], 'z')#'[100]')
            xyz_arrows.append(ar)
            ar = drawing.add_text([1*xyz_arrows_scale*100,-2.8*xyz_arrows_scale*100,-1.3*xyz_arrows_scale*100], 'y')#'[010]')
            xyz_arrows.append(ar)
            ar = drawing.add_text([0,0,2.9*xyz_arrows_scale*100], 'x')#'[001]')
            xyz_arrows.append(ar)

    if draw_stage:
        stage_dists = 2.2*65*np.array([18,20,20,22,25,28.25,30])-5
        stage_dists #*= scale
        stage_color = np.array([1.0,0.95,0.9,1])
        num_along_edge = 20
        stage_size = 35*100#*scale
        # stage rotator
        # top rotator
        ar = drawing.add_arrow([-stage_dists[0],0,0],
                               [-stage_dists[1],0,0],
                               head_fraction = 0, rod_radius = 0.4*stage_size/(stage_dists[1]-stage_dists[0]),
                               hat_radius = 0.1, nodes_per_circle = num_along_edge*4, color=stage_color, angle_offset=0)
        # handle dot
        ar = drawing.add_arrow([-0.5*(stage_dists[0]+stage_dists[1]),0.45*stage_size,0],
                               [-0.5*(stage_dists[0]+stage_dists[1]),0.4*stage_size,0],
                               head_fraction = 0, rod_radius = 17*0.3*(stage_dists[1]-stage_dists[0])/(0.55*stage_size),
                               hat_radius = 0.0001, nodes_per_circle = num_along_edge, color=stage_color, angle_offset=0)

    drawing.rot_x(phi)

    if draw_stage:
        # bottom rotator
        stage_color= np.copy(stage_color)
        stage_color[0:3]*= 0.8
        ar = drawing.add_arrow([-stage_dists[2],0,0],
                               [-stage_dists[3],0,0],
                               head_fraction = 0, rod_radius = 0.4*stage_size/(stage_dists[3]-stage_dists[2]),
                               hat_radius = 0.1, nodes_per_circle = num_along_edge*4, color=stage_color, angle_offset=0)
        # handle dot
        ar = drawing.add_arrow([-0.5*(stage_dists[2]+stage_dists[3]),0.45*stage_size,0],
                               [-0.5*(stage_dists[2]+stage_dists[3]),0.4*stage_size,0],
                               head_fraction = 0, rod_radius = 17*0.3*(stage_dists[1]-stage_dists[0])/(0.55*stage_size),
                               hat_radius = 0.0001, nodes_per_circle = num_along_edge, color=stage_color, angle_offset=0)

        nodes, faces = add_stage_section(stage_dists[4]+0.1, stage_dists[3], num_along_edge = num_along_edge, stage_size = stage_size)
        nodes[4*num_along_edge:,0] = -stage_dists[3]+0.1
        drawing.add_mesh(nodes, faces, facecolor=stage_color, abscolor=False)

        phi_ang = -phi*180/np.pi %360
        if phi_ang>180:
            phi_ang-=360
        if show_text:
            ar = drawing.add_text(np.copy(ar.node_locations[-1])+np.array([3,3,0]), f'𝜑 = {phi_ang:.2f}°')
        outstr.append(f'phi = {phi_ang:.2f} degrees')
    drawing.rot_z(chi)

    # mid part of stage
    if draw_stage:
        stage_color= np.copy(stage_color)
        stage_color[0:3]*= 1.2
        nodes, faces = add_stage_section(stage_dists[5]+0.1, stage_dists[4], num_along_edge = num_along_edge, stage_size = stage_size)
        temp = np.copy(nodes[:,1])
        nodes[:,1] = nodes[:,2]
        nodes[:,2] = temp
        drawing.add_mesh(nodes, faces, facecolor=stage_color, abscolor=False)
        chi_ang = -chi*180/np.pi %360
        if chi_ang>180:
            chi_ang-=360
        if show_text:
            ar = drawing.add_text(np.array([-stage_dists[4],0,-stage_size*0.5]), f'𝜒 = {-chi_ang:.2f}°')
        outstr.append(f'chi = {-chi_ang:.2f} degrees')


    drawing.rot_y(omega)
    # top part of stage
    if draw_stage:
        stage_color= np.copy(stage_color)
        stage_color[0:3]*= 0.8
        nodes, faces = add_stage_section(stage_dists[6], stage_dists[5], num_along_edge = num_along_edge, stage_size = stage_size)
        nodes[:4*num_along_edge,0] = -stage_dists[6]
        drawing.add_mesh(nodes, faces, facecolor=stage_color, abscolor=False)
        omega_ang = -omega*180/np.pi %360
        if omega_ang>180:
            omega_ang-=360
        if show_text:
            ar = drawing.add_text(np.array([-stage_dists[5],stage_size*0.5,0]), f'𝜔 = {-omega_ang:.2f}°')
        outstr.append(f'omega = {-omega_ang:.2f} degrees')

    # beam and lens
    k0 = apply_rotation(np.copy(params['k0_vector']), phi, chi, omega, beam_stage_angle)
    kh = apply_rotation(np.copy(params['kh_vector']), phi, chi, omega, beam_stage_angle)
    make_beam_and_lens(drawing, sample_node_collection.node_locations, kh, params['shape'],
            scale, extend_beam, params['beam_thickness'], params['transverse_width'],
            beam_opacity, draw_scattered_beam, extend_imaging_system, lens_scale)


    # translate arrows
    if draw_scattering_arrows:
        for ar in Q_arrows:
            ar*=30/8
            ar += np.array([-50,50,50])*q_arrows_scale*100
            ar.rotate_x(beam_stage_angle)

    if draw_axes:
        for ar in xyz_arrows:
            ar*=30/8
            ar += np.array([-10,-10,-10])*30/8*xyz_arrows_scale*100

    # rotate perspective nicely
    drawing.rot_y(-90*np.pi/180)
    drawing.rot_x(90*np.pi/180) # 95
    # move crystal_structure
    if not type(atom_list) == type(None):
        for atom in crystal_structure:
            atom += crysta_structure_position


    # draw
    drawing.draw_structure(ax)
    # enable render
    import plotly.graph_objects

    fig = plotly.graph_objects.Figure(data=ax['data'])
    camera = dict(
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=-0.2, y=0.4, z=2),
        projection = dict(
        type="orthographic",
        ),
    )
    scene = dict(
            xaxis = dict(visible=False, range=[-500000,500000]),
            yaxis = dict(visible=False, range=[-500000,500000]),
            zaxis = dict(visible=False, range=[-500000,500000]),
            #aspectmode='data',
            aspectratio = dict(x=80, y=80, z=80),
            annotations = ax['annotations'],

        )
    fig.update_layout( scene_camera = camera, scene = scene,
            margin=dict(t=0, r=0, l=0, b=0),
            height=800,
            width = 1100,
            )
    #for d in ax['data'][1:]:
    #    fig.add_trace(d)

    if not type(export_blender_filepath) == type(None):
        drawing.export_blender(export_blender_filepath)

    return drawing, fig, ax, outstr



class Drawing:
    '''
    maintains a list of nodes and objects
    Nodes are rotated
    Objects are transformed and rendered
    '''
    def __init__(self):
        self.objects = []
        self.nodes = []
    def add_text(self,location, text, color = [0,0,0,1], scale = 0.7):
        '''
        adds text
        input:
            location: coordinatex [x,y,z]
            text: the text
            color: text color
            scale: text scale
        '''
        text = object_classes.Text(location, text, color=color, scale=scale)
        self.nodes.append(text)
        self.objects.append(text)
        return text


    def add_mesh(self,node_locations, faces, meshcolor=[0,0,0,0], facecolor=[0,0.3,0.7,0.5], abscolor=False):
        '''
        adds a trianguar mesh, perhaps made by skimage.measure.marching_cubes or another member function, such as the make_arrow below
        input:
            node_loations: np.array of floats with size (n,3)
            faces: np.array of ints with size (m,3)
            meshcolor: color of the lines [r,g,b,a]
            facecolor: color of the faces [r,g,b,a]
        '''
        mesh = object_classes.Mesh(node_locations, faces, meshcolor=meshcolor, facecolor=facecolor, abscolor=abscolor)
        self.nodes.append(mesh)
        self.objects.append(mesh)
        return mesh

    def add_box(self, shape, center, meshcolor=[0,0,0,0], linewidth = 0.1, facecolor=[0,0.3,0.7,0.5], abscolor = False ):
        '''
        adds a box of shape at center
        input:
            shape: np array of length 3, shape in x, y z
            center:  np array of length 3, center
            meshcolor: float color of lines [r, g, b, a]
            linewidth: float linewidth
            facecolor: color of faces, [r, g, b, a]
            abscolor: bool
                True -> all faces will be facecolor
                False -> opacity and color of face modified according to face normal with respect to camera
        returns:
            nodes, faces, lines
        '''
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
        node_locations *= np.array(shape)
        node_locations += np.array(center)
        '''
        add tilt here
        '''
        # add nodes
        box_node_collection = object_classes.NodeCollection(node_locations)
        self.nodes.append(box_node_collection)
        # add faces
        if len(facecolor) < 4 or facecolor[3] > 0:
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                self.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = facecolor, abscolor=abscolor))
        # add lines
        if linewidth > 0:
            for seg in [[0,1], [1,2], [2,3], [3,0], [4,5], [5,6],[6,7], [7,4], [0,4], [1,5], [2,6], [3,7]]:
                self.objects.append(object_classes.BoxLine(box_node_collection, seg, meshcolor = meshcolor, linewidth = linewidth))
        return box_node_collection, self.objects[-(6+8):-8], self.objects[-8:]


    def add_atom_list(self, atom_list):
        '''
        adds atoms making a crystal structure
        input:
            atom_list = list of atoms or other objects
        returns:
            list of node objects (that support += and *= operators to move them in space)
        '''
        node_list = []
        for atom in atom_list:
            if hasattr(atom,'node_collection'):
                if not atom.node_collection in self.nodes:
                    self.nodes.append( atom.node_collection )
                    node_list.append( atom.node_collection )
            else:
                self.nodes.append( atom )
                node_list.append( atom )
            self.objects.append( atom )
        return node_list


    def draw_structure(self, ax):
        '''
        Renders the scene
        input:
            ax: plotly axes
        '''
        self.ax = ax
        # draw objects
        for obj in self.objects:
            obj.draw(self.ax)



    def rot_x(self,phi):
        '''
        rot_x: rotates around x axis (horisontal)
        input: phi, angle in radians
        '''
        for node in self.nodes:
            node.rotate_x(phi)

    def rot_y(self,phi):
        '''
        rot_y: rotates around y axis (vertical)
        input: phi, angle in radians
        '''
        for node in self.nodes:
            node.rotate_y(phi)

    def rot_z(self,phi):
        '''
        rot_z: rotates around z axis (out-of-plane)
        input: phi, angle in radians
        '''
        for node in self.nodes:
            node.rotate_z(phi)

    def add_arrow(self,start, end, head_fraction = 0.45, rod_radius = 0.1, hat_radius = 0.25, nodes_per_circle = 12, color=[0,0,1,1], angle_offset=0):
        '''
        adds an arrow based on a triangular mesh
        input:
            start: origin of arrow, array of lenght 3
            end: end of arrow, array of lenght 3
            head_fraction: optional float, fraction of arrow to be head
            rod_radius: optional float: radius of rod as fracton of length
            hat_radius: optional float: radius of hat as fracton of length
            nodes_per_circle: optional int: number of points along the circumfrence
            color: color of the faces [r,g,b,a]
            angle_offset: constant angle offset for the when nodes are placed in the list. Interacts with rendering in matplotlib to get correct face towards the camera
        '''
        mesh = make_arrow_mesh(start, end, head_fraction = head_fraction, rod_radius = rod_radius,
                               hat_radius = hat_radius, nodes_per_circle = nodes_per_circle, color=color, angle_offset=angle_offset)
        self.nodes.append(mesh)
        self.objects.append(mesh)

        return mesh

    def add_curved_arrow(self, origo, arrow_center, arrow_dir, angular_extent = np.pi/6, links = 5,
                         head_fraction = 0.25, rod_radius = 0.1, hat_radius = 0.15, nodes_per_circle = 12, color=[0,0,1,1], angle_offset=0):
        '''
        adds an arrow based on a triangular mesh
        input:
            origo: location [x,y,z], origin around which the arrow curves
            arrow_center: location [x,y,z],  center of the arrow, radius is defined by length of arrow_center-origo
            arrow_dir: vector [x,y,z]. The arrow curve in the plane np.cross(arrow_center-origo, arrow_dir), pointing in the direction where np.dot(arrow_dir, direction) is positive
            angular_extent: float, radians describing the extent for the arrow
            links: int, number of links (sylinders) for the arrow from base to head
            head_fraction: optional float, fraction of arrow to be head
            rod_radius: optional float: radius of rod as fracton of length
            hat_radius: optional float: radius of hat as fracton of length
            nodes_per_circle: optional int: number of points along the circumfrence
            color: color of the faces [r,g,b,a]
            angle_offset: constant angle offset for the when nodes are placed in the list. Interacts with rendering in matplotlib to get correct face towards the camera
        '''
        # parse
        origo = np. array(origo)
        arrow_center = np. array(arrow_center)
        radial_center = arrow_center-origo
        arrow_dir = np. array(arrow_dir)
        arrow_len = angular_extent*np.sqrt(np.sum((arrow_center-origo)**2)) # length along arrow
        rod_radius = rod_radius*arrow_len
        hat_radius = hat_radius*arrow_len
        plane_common_normal = np.cross(arrow_center-origo, arrow_dir)
        plane_common_normal = plane_common_normal/np.sqrt(np.sum(plane_common_normal**2))
        # declare node and facet vectors
        nodes = np.zeros((nodes_per_circle*(links+2)+2,3))
        facets = np.zeros((nodes_per_circle*(1+2*(links+1)+1),3), dtype=int)
        circles=[]
        # make butt of arrow
        radial = np.dot(rotation_matrix(plane_common_normal, -angular_extent/2), radial_center)
        nodes[0,:] = radial
        nodes[1:nodes_per_circle+1,:] = make_circ(plane_common_normal, normalize(radial), rod_radius, nodes_per_circle, angle_offset=angle_offset)+radial
        circles.append(np.arange(1,nodes_per_circle+1))
        facets[0:nodes_per_circle,:] = make_vertices_cone(0, circles[0] , nodes_per_circle)
        # make links
        angular_step = angular_extent*(1-head_fraction)/links
        for i in range(1, links+1):
            radial = np.dot(rotation_matrix(plane_common_normal, -angular_extent/2+angular_step*i), radial_center)
            #print(radial)
            circles.append(np.arange(1+nodes_per_circle*i, nodes_per_circle*(i+1)+1))
            nodes[i*nodes_per_circle+1:(i+1)*nodes_per_circle+1,:] = make_circ(plane_common_normal, normalize(radial), rod_radius, nodes_per_circle, angle_offset=angle_offset)+radial
            facets[(2*i-1)*nodes_per_circle:(2*i+1)*nodes_per_circle,:] = make_vertices_sylinder(circles[i-1], circles[i] , nodes_per_circle)
        # make hat brim
        i = links+1
        circles.append(np.arange(1+nodes_per_circle*i, nodes_per_circle*(i+1)+1))
        nodes[i*nodes_per_circle+1:(i+1)*nodes_per_circle+1,:] = make_circ(plane_common_normal, normalize(radial), hat_radius, nodes_per_circle, angle_offset=angle_offset)+radial
        facets[(2*i-1)*nodes_per_circle:(2*i+1)*nodes_per_circle,:] = make_vertices_sylinder(circles[i-1], circles[i] , nodes_per_circle)
        # make nose
        radial = np.dot(rotation_matrix(plane_common_normal, angular_extent/2), radial_center)
        nodes[-1,:] = radial
        facets[(2*i+1)*nodes_per_circle:,:] = make_vertices_cone(len(nodes)-1, circles[-1] , nodes_per_circle)
        # make mesh
        nodes[:,:]+=origo
        mesh = self.add_mesh(nodes, facets, facecolor=color)
        return mesh

    def export_blender(self, filepath):
        '''
        exports a list of all objects using pickle. Each object is a dictionary. This can be loaded into blender using an appropriate module
        iput:
            filepath: pickle file to store output
        returns:
            None
        '''
        top_level_list = []
        for obj in self.objects:
            top_level_list.append(obj.to_dict())
        import pickle
        pickle.dump(top_level_list, open( filepath, "wb" ))

    def add_lens(self, radius = 30, num_links = 4,
                r_curvature = 100, facing = np.array([1.0, 1.0, 1.0]),
                displacement = np.array([0,0,0])):

        def concave_surface(radius,num_links,r_curvature, facing, invert = False):
            branch_length = radius/num_links
            node_index = 0
            # 60 degree rotation matrix
            theta = np.pi/3
            rot_mat = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta),np.cos(theta)]])

            #start with center node
            nodes = [[np.array([0.0,0.0,0.0,node_index])]]
            facets = []
            # place firts ring
            nodes.append([])
            node_index += 1
            nodes[-1].append(nodes[-2][0]+np.array([0,-1,0,node_index]))

            dir = np.array([np.sqrt(3)/2,0.5,0,1])
            for i in range(5):

                temp_node = nodes[-1][-1]+dir
                node_index += 1
                temp_node[3] = node_index
                nodes[-1].append(temp_node)
                #print(dir[0:2],temp_node)

                facets.append(np.array([0,node_index-1,node_index]))
                dir[0:2] = np.matmul(rot_mat,dir[0:2])

            facets.append(np.array([0,1,node_index]))

            # we do num_links- circles out from center 7
            for i in range(num_links-1):
                nodes.append([]) # nodes stored as [[x,y,z,index]]
                # add node straight down
                node_index += 1
                nodes[-1].append(nodes[-2][0]+np.array([0,-1,0,0]))
                nodes[-1][-1][3] = node_index
                # start walking right and up
                dir = np.array([np.sqrt(3)/2,0.5,0,1])

                # walk around perrimiter
                l = len(nodes[-2])
                k = 0 #k walsk along inner nodes
                n = 0 #n counts turns
                while n<6:
                    # try step
                    temp_node = nodes[-1][-1]+dir
                    if all(np.abs( nodes[-2][(k+1)%l][0:2] - (nodes[-2][k%l][0:2]+dir[0:2]) ) < 0.01 ): #1.05*branch_length**2 > np.sum((temp_node[0:2]-nodes[-2][(k+1)%l][0:2])**2):
                        # i.e. new branch is parallel to the branch in the circle inside
                        node_index += 1
                        temp_node[3] = node_index
                        nodes[-1].append(temp_node)
                        facets.append(np.array([nodes[-2][k][3],node_index-1,node_index]))
                        facets.append(np.array([nodes[-2][(k+1)%l][3],nodes[-2][k][3],node_index]))
                        k+=1
                    else:
                        # we have moved past the edge
                        node_index += 1
                        temp_node[3] = node_index
                        nodes[-1].append(temp_node)
                        facets.append(np.array([nodes[-2][k%l][3],node_index-1,node_index]))
                        dir[0:2] = np.matmul(rot_mat,dir[0:2])
                        n+=1
            #unravel nodes
            nodes2 = []
            for i, lis in enumerate(nodes):
                if i==0:
                    nodes2.append(np.array([0,0,0]))
                    continue
                outer_nodes_start = len(nodes2)
                sqare_lim = radius/np.sqrt(2)
                for node in lis:
                    nodes2.append(node[0:3])
                    dist = num_links*branch_length
                    nodes2[-1] *= dist/np.sqrt(np.sum(node[0:3]**2))
                    if nodes2[-1][0] > sqare_lim: nodes2[-1][0]=sqare_lim
                    if nodes2[-1][0] < -sqare_lim: nodes2[-1][0]=-sqare_lim
                    if nodes2[-1][1] > sqare_lim: nodes2[-1][1]=sqare_lim
                    if nodes2[-1][1] < -sqare_lim: nodes2[-1][1]=-sqare_lim
                    nodes2[-1] *=i/num_links
                    nodes2[-1][2] = r_curvature-np.sqrt(r_curvature**2-np.sum(nodes2[-1][0:2]**2))
            outer_nodes = np.arange(outer_nodes_start,len(nodes2))
            nodes = np.array(nodes2)
            if invert:
                nodes[:,2]*=-1
            facets = np.array(facets, dtype = int)
            # rotate from facing a to facing b
            facing /= np.sqrt(np.sum(facing**2))
            original_facing = np.array([0, 0, 1.0])
            if not all(facing == original_facing):
                #https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
                v = np.cross(original_facing, [1.0,0,0])
                c = np.dot(original_facing, [1.0,0,0])
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                rotation_matix = np.eye(3) + vx + np.dot(vx,vx)/(1+c)
                nodes = np.matmul(rotation_matix, np.swapaxes(nodes,0,1))
                nodes = np.swapaxes(nodes,0,1)

                v = np.cross([1.0,0,0], facing)
                c = np.dot([1.0,0,0], facing)
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                rotation_matix = np.eye(3) + vx + np.dot(vx,vx)/(1+c)
                nodes = np.matmul(rotation_matix, np.swapaxes(nodes,0,1))
                nodes = np.swapaxes(nodes,0,1)
            return nodes, facets, outer_nodes

        nodes_0, facets_0, outer_nodes_0 = concave_surface(radius,num_links,r_curvature, facing) #front
        nodes_1, facets_1, outer_nodes_1 = concave_surface(radius,num_links,r_curvature, facing, invert = True) #back
        facets_1 += len(nodes_0)
        outer_nodes_1 += len(nodes_0)
        # make band connecting the two concave surfaces
        facets_2 = np.zeros((2*len(outer_nodes_0),3), dtype=int)
        l = len(outer_nodes_0)
        for i, _ in enumerate(outer_nodes_0):
            facets_2[2*i] = np.array([outer_nodes_0[i],outer_nodes_0[(i+1)%l],outer_nodes_1[i]])
            facets_2[2*i+1] = np.array([outer_nodes_1[i],outer_nodes_1[(i+1)%l],outer_nodes_0[(i+1)%l]])
        nodes = np.concatenate((nodes_0,nodes_1))
        facets = np.concatenate((facets_0,facets_1,facets_2))
        #for node in nodes:
            #rotate_x(node, phi)
            #rotate_z(node, chi)
            #rotate_y(node, omega)
            #rotate_y(node, -omega)
            #rotate_z(node, -chi)
            #rotate_x(node, -phi)
        nodes += np.array(displacement)
        mesh = object_classes.Mesh(nodes, facets, meshcolor=[0,0,0,0], facecolor=[0.9,0.9,0.9,1], abscolor=False)
        self.nodes.append(mesh)
        self.objects.append(mesh)
        return mesh

def normalize(vec):
    return vec/np.sqrt(np.sum(vec**2))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    thanks to user unutbu at stackexchange, https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
def make_circ(n2,n3,radius, nodes_per_circle, angle_offset=0):
    '''
    makes a circle of nodes in the n2 n3 plane
    '''
    circle_nodes = np.zeros((nodes_per_circle,3))
    for i in range(nodes_per_circle):
        ang = i*2*np.pi/nodes_per_circle + angle_offset
        circle_nodes[i] = radius*( n2*np.cos(ang) + n3*np.sin(ang))
    return circle_nodes
def make_vertices_cone(node0, circle_nodes, nodes_per_circle):
    '''
    makes a cone connecting node0 to nodes in circle_nodes
    '''
    facets = np.zeros((nodes_per_circle,3), dtype=int)
    for i in range(nodes_per_circle):
        facets[i] = [node0, circle_nodes[i], circle_nodes[(i+1)%nodes_per_circle]]
    return facets
def make_vertices_sylinder(circle_nodes_0, circle_nodes_1, nodes_per_circle):
    '''
    makes a sylinder connecting nodes in circle_nodes_0 to nodes in circle_nodes_1
    '''
    facets = np.zeros((2*nodes_per_circle,3), dtype=int)
    for i in range(nodes_per_circle):
        facets[2*i] = [circle_nodes_0[i], circle_nodes_1[i], circle_nodes_1[(i+1)%nodes_per_circle]]
        facets[2*i+1] = [circle_nodes_0[i], circle_nodes_0[(i+1)%nodes_per_circle], circle_nodes_1[(i+1)%nodes_per_circle]]
    return facets



def add_stage_section(dist_1, dist_2, num_along_edge = 30, stage_size = 12.5):
    '''
    Creates a triangular mesh describing a section of the stage, with curved top and bottom, described by dist_1 and dist_2
    The top curves along the second axis
    The bottom curves along the first axis
    input:
        dist_1: radial distance of inner circle
        dist_2: radial distance of outer circle
        num_along_edge: int, optional, number of nodes along each edge
        stage_size: size of the stage along a and b axis.
    '''
    edge = np.arange(-stage_size/2, stage_size/2+0.0001, stage_size/(num_along_edge-1))
    node_locations = np.zeros((num_along_edge*8,3))
    node_locations[0:4*num_along_edge,2] = np.sqrt(dist_1**2+stage_size**2/4)
    node_locations[4*num_along_edge:,2] = np.sqrt(dist_1**2+stage_size**2/4)
    # a coord
    node_locations[0:num_along_edge,1] = -stage_size/2
    node_locations[num_along_edge:2*num_along_edge,1] = edge
    node_locations[2*num_along_edge:3*num_along_edge,1] = stage_size/2
    node_locations[3*num_along_edge:4*num_along_edge,1] = -edge
    node_locations[4*num_along_edge:5*num_along_edge,1] = -stage_size/2
    node_locations[5*num_along_edge:6*num_along_edge,1] = edge
    node_locations[6*num_along_edge:7*num_along_edge,1] = stage_size/2
    node_locations[7*num_along_edge:,1] = -edge
    # b coord
    node_locations[0:num_along_edge,2] = edge
    node_locations[num_along_edge:2*num_along_edge,2] = stage_size/2
    node_locations[2*num_along_edge:3*num_along_edge,2] = -edge
    node_locations[3*num_along_edge:4*num_along_edge,2] = -stage_size/2
    node_locations[4*num_along_edge:5*num_along_edge,2] = edge
    node_locations[5*num_along_edge:6*num_along_edge,2] = stage_size/2
    node_locations[6*num_along_edge:7*num_along_edge,2] = -edge
    node_locations[7*num_along_edge:,2] = -stage_size/2
    # c coord
    node_locations[:4*num_along_edge,0] = -np.sqrt(dist_1**2-node_locations[:4*num_along_edge,1]**2)
    node_locations[4*num_along_edge:,0] = -np.sqrt(dist_2**2-node_locations[:4*num_along_edge,2]**2)

    faces = np.zeros((12*num_along_edge-4,3), dtype=int)
    '''# circumfrence
    for i in range(4*num_along_edge):
        faces.append([i,(i+1)%(4*num_along_edge),4*num_along_edge+i])
        faces.append([(i+1)%(4*num_along_edge),4*num_along_edge+(i+1)%(4*num_along_edge),4*num_along_edge+i])'''
    i_edge = np.arange(4*num_along_edge, dtype=int)
    faces[0:4*num_along_edge,0] = i_edge
    faces[0:4*num_along_edge,1] = (i_edge+1)%(4*num_along_edge)
    faces[0:4*num_along_edge,2] = 4*num_along_edge+i_edge
    faces[4*num_along_edge:8*num_along_edge,0] = (i_edge+1)%(4*num_along_edge)
    faces[4*num_along_edge:8*num_along_edge,1] = 4*num_along_edge+(i_edge+1)%(4*num_along_edge)
    faces[4*num_along_edge:8*num_along_edge,2] = 4*num_along_edge+i_edge

    '''# top
    for i in range(num_along_edge,2*num_along_edge-1):
        faces.append([i,i+1,5*num_along_edge-i-1])
        faces.append([i+1,5*num_along_edge-i-2,5*num_along_edge-i-1])
    '''
    i_edge = np.arange(num_along_edge,2*num_along_edge-1, dtype=int)
    faces[8*num_along_edge:9*num_along_edge-1,0] = i_edge
    faces[8*num_along_edge:9*num_along_edge-1,1] = i_edge+1
    faces[8*num_along_edge:9*num_along_edge-1,2] = 5*num_along_edge-i_edge-1
    faces[9*num_along_edge-1:10*num_along_edge-2,0] = i_edge+1
    faces[9*num_along_edge-1:10*num_along_edge-2,1] = 5*num_along_edge-i_edge-2
    faces[9*num_along_edge-1:10*num_along_edge-2,2] = 5*num_along_edge-i_edge-1
    '''
    # bottom
    for i in range(num_along_edge-1):
        faces.append([i+4*num_along_edge,i+1+4*num_along_edge,7*num_along_edge-i-1])
        faces.append([i+1+4*num_along_edge,7*num_along_edge-i-2,7*num_along_edge-i-1])
    faces = np.array(faces)'''
    i_edge = np.arange(num_along_edge-1, dtype=int)
    faces[10*num_along_edge-2:11*num_along_edge-3,0] = i_edge+4*num_along_edge
    faces[10*num_along_edge-2:11*num_along_edge-3,1] = i_edge+1+4*num_along_edge
    faces[10*num_along_edge-2:11*num_along_edge-3,2] = 7*num_along_edge-i_edge-1
    faces[11*num_along_edge-3:12*num_along_edge-4,0] = i_edge+1+4*num_along_edge
    faces[11*num_along_edge-3:12*num_along_edge-4,1] = 7*num_along_edge-i_edge-2
    faces[11*num_along_edge-3:12*num_along_edge-4,2] = 7*num_along_edge-i_edge-1

    return node_locations, faces

def make_arrow_mesh(start, end, head_fraction = 0.45, rod_radius = 0.1, hat_radius = 0.25, nodes_per_circle = 12, color=[0,0,1,1], angle_offset=0):
    '''
    creates a triangular mesh in the shape of an arrow
    input:
        start: origin of arrow, array of lenght 3
        end: end of arrow, array of lenght 3
        head_fraction: optional float, fraction of arrow to be head
        rod_radius: optional float: radius of rod as fracton of length
        hat_radius: optional float: radius of hat as fracton of length
        nodes_per_circle: optional int: number of points along the circumfrence
        color: color of the faces [r,g,b,a]
        angle_offset: constant angle offset for the when nodes are placed in the list. Interacts with rendering in matplotlib to get correct face towards the camera
    '''
    # create basis of orthogonal normal vectors, n1 along arrow
    start = np. array(start)
    end = np. array(end)
    v1 =  end-start # along arrow
    n1 = v1/np.sqrt(np.sum(v1**2))
    v2 = np.cross(n1,np.array([1,0,0]))
    if np.sum(np.abs(v2)) < np.sum(np.abs(np.cross(n1,np.array([0,1,0])))):  # if v1 is parallel to [1,0,0]
        v2 = np.cross(n1,np.array([0,1,0]))
    n2 = v2/np.sqrt(np.sum(v2**2))
    n3 = np.cross(n1,n2)
    nodes = np.zeros((nodes_per_circle*3+2,3))
    d_angle = -2*np.pi/nodes_per_circle/2
    nodes[1:nodes_per_circle+1,:] = make_circ(n2,n3,rod_radius, nodes_per_circle, angle_offset=angle_offset)
    nodes[nodes_per_circle+1:2*nodes_per_circle+1,:] = make_circ(n2,n3,rod_radius, nodes_per_circle, angle_offset=angle_offset+d_angle) + (1-head_fraction)*n1
    nodes[2*nodes_per_circle+1:3*nodes_per_circle+1,:] = make_circ(n2,n3,hat_radius, nodes_per_circle, angle_offset=angle_offset+2*d_angle) + (1-head_fraction)*n1
    nodes[3*nodes_per_circle+1,:] = n1
    nodes = nodes * np.sqrt(np.sum(v1**2))
    nodes += start

    facets = np.zeros((nodes_per_circle*(1+2+2+1),3), dtype=int)
    circ_0 = np.arange(1,nodes_per_circle+1)
    circ_1 = np.arange(nodes_per_circle+1,2*nodes_per_circle+1)
    circ_2 = np.arange(2*nodes_per_circle+1,3*nodes_per_circle+1)
    facets[0:nodes_per_circle,:] = make_vertices_cone(0, circ_0 , nodes_per_circle)
    facets[nodes_per_circle:3*nodes_per_circle,:] = make_vertices_sylinder(circ_0, circ_1, nodes_per_circle)
    facets[3*nodes_per_circle:5*nodes_per_circle,:] = make_vertices_sylinder(circ_1, circ_2, nodes_per_circle)
    facets[5*nodes_per_circle:6*nodes_per_circle,:] = make_vertices_cone(3*nodes_per_circle+1, circ_2, nodes_per_circle)

    return object_classes.Mesh(nodes, facets, meshcolor=[0,0,0,0], facecolor=color, abscolor=False)


def vector_to_hkl(vector):
        '''
        Converts a vector, i.e. [0.5,0.33,0.33] to hkl, i.e. [3,2,2]
        input:
            vector: np.array of length 3
        '''
        vector2 = vector/np.max(vector)
        vector2[vector2<10**-4] = 0
        # we assume the H, K, and L are factors of 2*2*2*3*3*5*7*11*13
        vector2 *= 2*2*2*3*3*5*7*13*11
        # at this point we assume vector2 only contains ints
        int_vector = np.array(vector2+0.5, dtype='int') # +0.5 for correct rounding for pos numbers
        int_vector[int_vector<0] -= 1 # correct rounding for neg numbers
        int_vector = int_vector//np.gcd.reduce(np.abs(int_vector)) # divide by the greatest common divisor
        return int_vector



def rotate(x,z,phi):
    '''
    common rotation function, rotates
    rotation in the plane x,z according to phi
    input:
        x,z floats or numpy arrays
        input: phi, angle in radians
    returns:
        rotated x, rotated z
    '''
    return x*np.cos(phi)-z*np.sin(phi),z*np.cos(phi)+x*np.sin(phi)

def rotate_x(loc,phi):
    '''
    rotate_x: rotates around x axis (horisontal)
    input: phi, angle in radians
    '''
    x,z=rotate(loc[1],loc[2],phi)
    loc[1] = x
    loc[2] = z
def rotate_y(loc,phi):
    '''
    rotate_y: rotates around y axis (vertical)
    input: phi, angle in radians
    '''
    x,z=rotate(loc[0],loc[2],phi)
    loc[0] = x
    loc[2] = z
def rotate_z(loc,phi):
    '''
    rotate_z: rotates around z axis (out-of-plane)
    input: phi, angle in radians
    '''
    x,z=rotate(loc[0],loc[1],phi)
    loc[0] = x
    loc[1] = z

def get_phi_chi_omega(k0, kh, tilt_to_vector ):

    phi = np.arctan2(k0[1],k0[2])
    #print(tilt_to_vector)
    vec = np.copy(tilt_to_vector)
    rotate_x(vec,phi)
    chi = np.arctan2(vec[1],np.sqrt(vec[0]**2+vec[2]**2))
    rotate_z(vec,chi)
    omega = -np.arctan2(vec[2],vec[0])+np.pi
    vec2 = apply_rotation(tilt_to_vector, phi, chi, omega, 0)
    beam_stage_angle = 0
    def rot_to_fit(inp):
        phi = inp[0]
        chi = inp[1]
        omega = inp[2]
        beam_stage_angle = inp[3]
        tt_rot = normalize(apply_rotation(np.copy(tilt_to_vector), *inp))
        k0_rot = normalize(apply_rotation(np.copy(k0), *inp))
        kh_rot = normalize(apply_rotation(np.copy(kh), *inp))
        #print(vec, vec2)
        #print( np.abs(vec[0]-1)*1000+np.abs(vec[2])+np.abs(vec[1])+np.abs(vec2[0])+np.abs(vec2[1]))
        k0_e = np.abs(k0_rot[0])+np.abs(k0_rot[1])+k0_rot[2]
        kh_e = np.abs(kh_rot[1])
        tt_e = np.abs(tt_rot[2]-1)
        return k0_e + kh_e + tt_e
        return np.abs(vec[0]-1)+np.abs(vec2[2]-1)+np.abs(vec[2])+np.abs(vec[1])+np.abs(vec2[0])+np.abs(vec2[1])# + np.abs(params['k0_vector'][0]) + np.abs(params['k0_vector'][1])

    out= scipy.optimize.minimize(rot_to_fit,np.array([phi,chi,omega, beam_stage_angle]), tol=10**-10)
    #print(out)
    phi = out['x'][0]
    chi = out['x'][1]
    omega = out['x'][2]
    beam_stage_angle = out['x'][3]
    return phi, chi, omega, beam_stage_angle

def apply_rotation(vec, phi, chi, omega, beam_stage_angle):
    rotate_x(vec, phi)
    rotate_z(vec, chi)
    rotate_y(vec, omega)
    rotate_x(vec, beam_stage_angle)
    return vec


def make_beam_and_lens(drawing, sample_nodes, kh, shape, scale, extend_beam, beam_thickness, transverse_width, beam_opacity, draw_scattered_beam, extend_imaging_system, lens_scale):
    # add beam
    if not type(extend_beam) == type(None):
        node_locs = np.zeros((8,3))
        node_locs[[0,1,4,5],0] += 0.5*beam_thickness*scale
        node_locs[[2,3,6,7],0] -= 0.5*beam_thickness*scale
        node_locs[[0,3,4,7],1] += 0.5*transverse_width*scale
        node_locs[[2,1,6,5],1] -= 0.5*transverse_width*scale
        node_locs[0:4,2] -= extend_beam[0]*scale
        node_locs[4:,2]  += extend_beam[1]*scale
        beam_node_locs = np.copy(node_locs)
        box_node_collection = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(box_node_collection)
        if 1:
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,beam_opacity], abscolor=False))
        #drawing.objects.append(object_classes.BoxFacet(box_node_collection, [0,1,2,3], facecolor = [1,0.6,0,beam_opacity*2], abscolor=False))
        #drawing.objects.append(object_classes.BoxFacet(box_node_collection, [0,1,2,3], facecolor = [1,0.6,0,beam_opacity*2], abscolor=False))
    sample_norm = np.cross(sample_nodes[1]-sample_nodes[0],sample_nodes[2]-sample_nodes[0])
    lens_center = normalize(kh)*extend_imaging_system[1]
    # scattered beam to lens
    if draw_scattered_beam:
        node_locs = np.zeros((8,3))
        node_locs[0] = line_plane_collision(sample_norm, sample_nodes[0], beam_node_locs[0], beam_node_locs[4])
        node_locs[1] = line_plane_collision(sample_norm, sample_nodes[0], beam_node_locs[1], beam_node_locs[5])
        node_locs[2] = line_plane_collision(sample_norm, sample_nodes[4], beam_node_locs[2], beam_node_locs[6])
        node_locs[3] = line_plane_collision(sample_norm, sample_nodes[4], beam_node_locs[3], beam_node_locs[7])
        for i in range(4,8):
            node_locs[i] = line_plane_collision(kh, lens_center, node_locs[i-4], node_locs[i-4]+kh)
        sample_in_lens = np.copy(node_locs[4:])
        # make object
        box_node_collection = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(box_node_collection)
        if 1: # projection along kh
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,0.25], abscolor=False))
        node_locs[4:] = 0.01*node_locs[4:]+0.99*np.average(node_locs[4:],axis = 0)
        box_node_collection = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(box_node_collection)
        if lens_scale>0:
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,beam_opacity], abscolor=False))
    # lens
    if lens_scale>0:
        lens_radius = lens_scale*0.5*np.sqrt(2)*scale*5000#np.max([params['shape'][0],params['shape'][1]])
        r_curvature = lens_radius*3
        delta_lens = normalize(kh)
        delta_lens*= 1.1*2*(r_curvature-np.sqrt(r_curvature**2-lens_radius**2))
        num_lenses = 10
        for i in np.arange(num_lenses):
            i-= (num_lenses-1)/2-1
            mesh = drawing.add_lens(radius = lens_radius, num_links = 15, r_curvature = r_curvature,
                facing = -kh, displacement = lens_center + delta_lens*i)
    # beam behind lens
    if lens_scale>0:
        magnification = 3
        node_locs = np.zeros((8,3))
        node_locs[:4] = sample_in_lens
        node_locs[4:] = sample_in_lens*3+lens_center
        # make object
        box_node_collection = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(box_node_collection)
        if 1: # projection along kh
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,0.25], abscolor=False))
        node_locs[:4] = 0.01*node_locs[:4]+0.99*np.average(node_locs[:4],axis = 0)
        box_node_collection = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(box_node_collection)
        if lens_scale>0:
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,beam_opacity], abscolor=False))

    return


def line_plane_collision(plane_norm, plane_point, line_start, line_end):
    direction = line_end-line_start
    ndotu = plane_norm.dot(direction)
    w = line_start - plane_point
    si = -plane_norm.dot(w) / ndotu
    Psi = w + si * direction + plane_point
    return Psi
