import numpy as np
import copy
import three_d_draw_object_classes as object_classes
import time
import streamlit as st

def set_backend(backend):
    '''
    sets the backend in object_classes
    input:
        backend: string 'matplotlib' or 'mayavi'
    '''
    object_classes.set_backend(backend)
    return

t0=time.time()
timing=False
timing_render=False
def print_time(text=''):
    if timing:
        global t0
        print(round((time.time()-t0)*100)/100, text)
        t0 = time.time()
def print_time_render(text=''):
    if timing_render:
        global t0
        print(round((time.time()-t0)*100)/100, text)
        t0 = time.time()
def set_edge_to_const(mat, border_val):
    '''
    sets the outer edges of mat to value
    input:
        mat: 3d numpy array
        border_val: float
    returns:
        mat
    '''
    mat[0,:,:] = border_val
    mat[:,0,:] = border_val
    mat[:,:,0] = border_val
    mat[-1,:,:] = border_val
    mat[:,-1,:] = border_val
    mat[:,:,-1] = border_val
    return mat

def make_3d_perspective(ex,
                fig=None, ax=None, export_blender_filepath=None,
                # deals with rendering the figure
                view_angle = 3, fig_size=None, tilt_to_vector = None,
                factor = 2000, final_rotation_function = None,
                scale = 0.1,
                # deals with additional objects to render
                volume_data = None,
                atom_list = None, crysta_structure_position = [-70, -10, 60],
                beam_draw_mesh = True,
                extend_beam = (1,0), beam_step_size = 1, beam_opacity = 1,
                extend_imaging_system=(0,1), imaging_system_lw = 1,
                draw_scattered_beam = True, scattered_beam_lw =0,
                lens_scale = 0,
                bounding_box_facecolor = [0,0,0,0.3],
                # image at the end of the imaging system
                intensity_imaging_system = None,
                intensity_imaging_system_opacity = None,
                # detector
                detector_map = None,
                #extend_detector = None,
                #detector_intensity = None,
                show_beam_at_end = False,
                # stages and arrows
                draw_stage = False,
                draw_axes = True,
                draw_scattering_arrows = True,
                draw_face_axes = False, crystal_rotation_function = None,
                draw_curved_arrows = False,
                # misc
                img_colormap = 'viridis',
                arrow_nodes_per_circle=12,
                ):
    '''
    Makes a matplotlib figure illustrating the experiment
    input:
        ex: DFXRM experiment object
        fig: optional, matplotlib or mayavi figure. A new figure is created if None
        ax: matlotlib axis, or None
        scale: float default = 1.0, multiplication factor for the simulated volume
        volume_data: optional list of meshes to be visualized. list of object_classes.Mesh objects
        beam_step_size: optional int,  step size for calculating triangle mesh of the beam
        beam_opacity: opacity of the beam
        intensity_imaging_system: 2d image data to show at the end of the the projected imaging system (intensity corresponding to ex.Eh_detector)
        intensity_imaging_system_opacity: opacity data of the above used in blender
        show_beam_at_end: optional bool, if True, show a pcolormesh of beam at end'
        view_angle: mayavi view angle, optional
        fig_size: size of the figure, in inches for matplotlib, or pixels for mayavi
        tilt_to_vector: vector of lenght 3 that defines 'up' or none
        extend_beam: tuple or None, default (1,0). The beam is extended in (before, after) directions along the [001] axis equal to this value times the cell lenght along [001]
                    If None, the beam is not shown
        extend_imaging_system: tuple (forward,back), default = (5,5). The integrated field projection is extended in both directions along the [001] axis equal to this value times the cell lenght along [001].
                        It is invisible ulness imaging_system_lw is also set.
                        also determines the posistion of lenses, intensity_imaging_system, and detector_map
        imaging_system_lw: float, default is 0, linewith of the box of extended intensity
        draw_scattered_beam: bool default True, if true draws the scattered beam
        scattered_beam_lw: linewidth to draw the edges of the scattered beam
        extend_detector: float or None, default = None. The detector projection is extended in both directions along the [001] axis equal to this value times the cell lenght along [001]
        export_blender_filepath: string, optional filepath for blender export. Default None
        draw_stage: bool default False. If True, draw the stage
        draw_axes: bool default True. If False, do not draw the arrows indicating the sample geometry
        draw_scattering_arrows: bool default True. If False, do not draw the arrows indicating the scattering vector
        draw_face_axes: bool default False. If True, draw the arrows indicating the sample geometry at the facets
        draw_curved_arrows: bool default True. If False, do not draw the arrows corresponding to rotation
        arrow_nodes_per_circle: int default 12. Nodes per circel when drawing arrwos
        atom_list: None or list of atoms and other objects to render crystal structure
        crysta_structure_position: position to render the crystal strcuture (with respect to the volume center of the simulated volume. This is in the coordinate system of the drawing: [horizontal, vertical, out-of-plane]
        factor: float, default 200. Distance from camera to structure. High number: parallel projection
        final_rotation_function: optional default None. If not None, should be a function like:
            def final_rotation_function(drawing):
                drawing.rot_y(-90*np.pi/180) # y in coordinate system of drawing is vertical
                drawing.rot_x(100*np.pi/180) # x in coordinate system of drawing is horizontal
                drawing.rot_z(-20*np.pi/180) # z in coordinate system of drawing is out of plane
                return
        img_colormap: string, default 'viridis', colormap to use when rendering images, a sting as accepted by mayavi/matplotlib
        lens_scale: float, default=0, scale to render the lens in. If 0, do not show the lens
        detector_map: map containing 'pixel_size', 'y_axis', 'x_axis' and 'intensity' draws the image on the detector at the location given by extend_imaging_system
        crystal_rotation_function: function tha takes an array of lenght three describing a vector in the crystal basis, and returns the vector in simulation coordiantes
    returns:
        drawing, fig, ax
    '''

    if timing:
        global t0
        t0 = time.time()
    if object_classes.backend == 'matplotlib':
        print(' Using matplotlib backend, objects are rendered as layers, and do not intersect correctly.')
        if type(fig) == type(None):
            if type(fig_size) == type(None):
                fig_size = (12,10)
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1,1,figsize=fig_size, dpi=100)
    elif object_classes.backend == 'mayavi':
        if type(fig) == type(None):
                import mayavi.mlab
                mayavi.mlab.options.offscreen = True
                if type(fig_size) == type(None):
                    fig_size = (12,6)
                fig_size = (fig_size[0]*200, fig_size[1]*200)
                fig = mayavi.mlab.figure(bgcolor=(1,1,1), fgcolor=None, engine=None, size=fig_size)
                fig.scene.disable_render = True # supress rendering after each object is added. Set to False at end
                ax = fig
    elif object_classes.backend == 'plotly':
        ax = {} # <- this list will be populated by plotly objects
        ax['data'] = []
        ax['annotations'] = []
        fig = None
    print_time('# add drawing')
    drawing = Drawing()
    # add crystal structure if cif specified
    if not type(atom_list) == type(None):
        print_time('# add crystal structure')
        crystal_structure = drawing.add_atom_list(atom_list)
    # add simulated box
    print_time('# add simulated box')
    box_shape = ex.shape*ex.del_x*scale
    box_nodes, _, __ = drawing.add_box(box_shape, [0,0,0], meshcolor=[0,0,0,1], linewidth = 2, facecolor=bounding_box_facecolor, abscolor = False)
    # text box size
    print_time('# text box size')
    drawing.add_text([0,1+0.5*box_shape[1],1+0.5*box_shape[2]], f'{box_shape[0]/scale:.2f} um')
    drawing.add_text([1+0.5*box_shape[0],0,1+0.5*box_shape[2]], f'{box_shape[1]/scale:.2f} um')
    drawing.add_text([-3-0.5*box_shape[0],0.5*box_shape[1],-8], f'{box_shape[2]/scale:.2f} um')

    # beam at end
    print_time('# beam at end')
    if show_beam_at_end:
        beam_at_end = np.abs( np.fft.ifft2( np.fft.fft2(ex.E0_entrance[...,0]) * np.exp(-2j*np.pi*ex.L/ex.k0_vector[2]*(ex.qx*ex.k0_vector[0] + ex.qy*ex.k0_vector[1])) )  )**2
        pcolormesh = object_classes.QuadrilateralColormesh(box_nodes, [4,5,6,7], beam_at_end, colormap = img_colormap)
        drawing.objects.append(pcolormesh)
    # add projection of integration
    #if not type(extend_imaging_system) == type(None):
    print_time('# add projection of integrated field')
    if 1: # always make this, but sometimes it is hidden depending on imaging_system_lw
        det_projection_nodes, _, __ = drawing.add_box(box_shape*np.array([1,1,0]), [0,0,0],
                                                        meshcolor=[0,0,0,1], linewidth = imaging_system_lw, facecolor=[0,0,0,0], abscolor = False)
        node_locs = det_projection_nodes.node_locations
        kh_norm = ex.kh_vector/np.sqrt(np.sum(ex.kh_vector**2))
        node_locs[:,:] -= np.sum(node_locs*kh_norm[np.newaxis,:], axis=1)[:,np.newaxis]*kh_norm[np.newaxis,:]
        node_locs[:4,:] -= box_shape[2]*(1+extend_imaging_system[0])*kh_norm[np.newaxis,:]/kh_norm[2]
        node_locs[4:,:] += box_shape[2]*extend_imaging_system[1]*kh_norm[np.newaxis,:]/kh_norm[2]
        node_locs[:,2] += 0.5*box_shape[2]
    if draw_scattered_beam:
        # add projection from beam intersection
        #z_vec = np.array([0,0,1])
        #beam_y = np.cross(ex.Q_vector,z_vec)
        #beam_y = beam_y/np.sqrt(np.sum(beam_y**2))
        node_locs = np.zeros((8,3))

        surface_angle = -np.arctan2(ex.k0_vector[1], -ex.k0_vector[0])
        mid_point = ex.params['BEAM']['mid'] # mid in relative units
        mid_point = [mid_point[0] * ex.shape[0]*ex.del_x[0], mid_point[1] * ex.shape[1]*ex.del_x[1]]

        # front intesection of beam
        beam_loc_0 = np.array([0,1,0])*ex.params['BEAM']['transverse_width']
        node_locs[0,0] = np.cos(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]+beam_loc_0[0] +mid_point[0]) - np.sin(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]+beam_loc_0[1] +mid_point[1])
        node_locs[0,1] = np.cos(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]+beam_loc_0[1] +mid_point[1]) + np.sin(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]+beam_loc_0[0] +mid_point[0])
        node_locs[0,:] *= scale
        node_locs[1,0] = np.cos(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]-beam_loc_0[0] +mid_point[0]) - np.sin(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]-beam_loc_0[1] +mid_point[1])
        node_locs[1,1] = np.cos(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]-beam_loc_0[1] +mid_point[1]) + np.sin(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]-beam_loc_0[0] +mid_point[0])
        node_locs[1,:] *= scale

        node_locs[1,:] = node_locs[1,:]
        # back intersection of beam
        node_locs[2,:] = node_locs[1,:] + box_shape[2]*ex.k0_vector/ex.k0_vector[2]
        node_locs[3,:] = node_locs[0,:] + box_shape[2]*ex.k0_vector/ex.k0_vector[2]
        # upper projected nodes
        node_locs[4,:] = node_locs[0,:] + box_shape[2]*ex.kh_vector/ex.kh_vector[2]
        node_locs[5,:] = node_locs[1,:] + box_shape[2]*ex.kh_vector/ex.kh_vector[2]
        # bottom projected nodes
        node_locs[6,:] = node_locs[1,:] + box_shape[2]*ex.k0_vector/ex.k0_vector[2]
        node_locs[7,:] = node_locs[0,:] + box_shape[2]*ex.k0_vector/ex.k0_vector[2]
        # extend back end of box to the detector
        node_locs[4:,:] += box_shape[2]*extend_imaging_system[1]*ex.kh_vector[np.newaxis,:]/ex.kh_vector[2]
        node_locs[:,:] += -0.5*box_shape[2]*np.array([0,0,1])
        # make object
        box_node_collection = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(box_node_collection)
        # add faces
        if 1:
            for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,beam_opacity], abscolor=False))
        # add lines
        if scattered_beam_lw >0:
            for seg in [[0,1], [1,2], [2,3], [3,0], [4,5], [5,6],[6,7], [7,4], [0,4], [1,5], [2,6], [3,7]]:
                drawing.objects.append(object_classes.BoxLine(box_node_collection, seg, meshcolor = [0,0,0,1], linewidth = scattered_beam_lw))

    if not type(detector_map) == type(None):
        print_time('# add detector')
        extend_detector = extend_imaging_system[1]
        node_locs = np.array(
            [[0, 0, 0.5*box_shape[2]*(1+2*extend_detector)]]*4,
            )
        #print(node_locs)

        # project along kh
        node_locs[:,1] +=  (node_locs[:,2]-0.5*box_shape[2])*ex.kh_vector[1]/ex.kh_vector[2] # *np.cos((ex.theta-angle_Q_z_axis))*np.sin(angle_Q_projceted_XY_wrtY)
        node_locs[:,0] +=  (node_locs[:,2]-0.5*box_shape[2])*ex.kh_vector[0]/ex.kh_vector[2]# *np.cos((ex.theta-angle_Q_z_axis))*np.cos(angle_Q_projceted_XY_wrtY)
        node_locs /= scale
        # move nodes according to x and y axis of detector
        detector_shape = np.array([detector_map['x_axis']*detector_map['pixel_size'][0]*detector_map['intensity'].shape[0],
                                   detector_map['y_axis']*detector_map['pixel_size'][1]*detector_map['intensity'].shape[0]])
        if 'scale' in detector_map.keys():
           detector_shape *= detector_map['scale']
        if 'render_scale' in detector_map:
            detector_shape = detector_map['render_scale']*detector_shape
        node_locs[0] += 0.5*(-detector_shape[0]-detector_shape[1])
        node_locs[1] += 0.5*( detector_shape[0]-detector_shape[1])
        node_locs[2] += 0.5*( detector_shape[0]+detector_shape[1])
        node_locs[3] += 0.5*(-detector_shape[0]+detector_shape[1])
        node_locs*= scale
        node_locs = object_classes.NodeCollection(node_locs)
        drawing.nodes.append(node_locs)
        pcolormesh = object_classes.QuadrilateralColormesh(node_locs, [0,1,2,3],
                detector_map['intensity'], opacity_data = intensity_imaging_system_opacity, colormap = img_colormap)
        drawing.objects.append(pcolormesh)

    '''# add projection of detector
    if not type(extend_detector) == type(None):
        ## get detector geometry
        print_time('# add projection of detector')
        fov_center = np.array(ex.params['DETECTOR']['fov_center'])
        shape_detec = np.array(ex.params['DETECTOR']['shape'], dtype='int')
        M = float(ex.params['DETECTOR']['magnif'])
        del_x_detec = np.array(ex.params['DETECTOR']['pix_size']) / M
        detector_shape = (shape_detec[0]*del_x_detec[0], shape_detec[1]*del_x_detec[1])
        invm = np.linalg.inv(ex.M_imag)
        # make detector box
        det_detector_nodes, _, __ = drawing.add_box((1,1,1), (0,0,0), angle_vector=np.array([0,0,0]), meshcolor=[0,0,0,1], linewidth = 0.2, facecolor=[0,0,0,0], abscolor = False)
        det_detector_nodes.node_locations *= scale
        node_locs = det_detector_nodes.node_locations
        # move nodes in box in the x-y plane
        det_x_r = invm[0,0] * node_locs[:,0] + invm[0,1] * node_locs[:,1] # r = rotated and transformed
        det_y_r = invm[1,0] * node_locs[:,0] + invm[1,1] * node_locs[:,1]
        node_locs[:,0] = (det_x_r*detector_shape[0]+ex.shape[0]*(fov_center[0]-0.5))*ex.del_x[0] # rescale and center detector projection
        node_locs[:,1] = (det_y_r*detector_shape[1]+ex.shape[1]*(fov_center[1]-0.5))*ex.del_x[1]
        node_locs[:,2] = node_locs[:,2] * ex.shape[2] * ex.del_x[2] * (1+2*extend_detector)
        # transpose according to angle of kh_vector
        print_time('# transpose according to angle of kh_vector')
        node_locs[:,1] +=  (node_locs[:,2]-0.5*box_shape[2])*ex.kh_vector[1]/ex.kh_vector[2]
        node_locs[:,0] +=  (node_locs[:,2]-0.5*box_shape[2])*ex.kh_vector[0]/ex.kh_vector[2]
        # move nodes to the plane of the detector (i.e the ex.Q_vector-np.cross(ex.kh_vector, ex.Q_vector ) plane)
        detector_normal = np.cross( ex.Q_vector, np.cross( ex.kh_vector, ex.Q_vector ))
        detector_normal = detector_normal/np.sqrt(np.sum(detector_normal**2))
        k_h_norm = ex.kh_vector/np.sqrt(np.sum(ex.kh_vector**2))
        for plane in [0,4]:
            for i in [0,2,3]:
                dist = np.dot((node_locs[i+plane,:]-node_locs[1+plane,:]),detector_normal)
                node_locs[i+plane,:] -= dist*k_h_norm'''
    # add beam
    if not type(extend_beam) == type(None):
        print_time('# add beam')
        if beam_draw_mesh:
            input_beam = np.zeros((*(ex.E0_entrance[:,:,0].shape),beam_step_size*2))
            for i in range(input_beam.shape[-1]):
                input_beam[:,:,i]= ex.E0_entrance[:,:,0].real
            beam_level = 0.5*np.max(input_beam)
            # beam is splitt to improve rendering for the matplotlib backend, but is not needed for mayavi, but do it anyway ¯\_(ツ)_/¯
            # beam outside (enterence)
            mesh  = drawing.add_isosurface_mesh(isosurface_mesh_from_matrix(input_beam, level = beam_level,
                                    step_size = beam_step_size, del_x = ex.del_x, meshcolor=[0,0,0,0], facecolor=[1,0.6,0,beam_opacity], abscolor=True))
            vertices = mesh.node_locations
            vertices[:,2] -= np.average([np.max(vertices[:,2]),np.min(vertices[:,2])])
            vertices[:,2] *= 0.5*ex.del_x[2]*ex.shape[2]/np.max(vertices[:,2])
            vertices[:,2] =  (vertices[:,2]*extend_beam[0]-(0.5+extend_beam[0]/2)*ex.del_x[2]*ex.shape[2]) # shift to front
            vertices[:,1] += (vertices[:,2] + 0.5*ex.del_x[2]*ex.shape[2])*ex.k0_vector[1]/ex.k0_vector[2]
            vertices[:,0] += (vertices[:,2] + 0.5*ex.del_x[2]*ex.shape[2])*ex.k0_vector[0]/ex.k0_vector[2]
            vertices *= scale
            # beam inside
            mesh  = drawing.add_isosurface_mesh(isosurface_mesh_from_matrix(input_beam, level = beam_level,
                                    step_size = beam_step_size, del_x = ex.del_x, meshcolor=[0,0,0,0], facecolor=[1,0.6,0,beam_opacity], abscolor=True))
            vertices = mesh.node_locations
            vertices[:,2] -= np.average([np.max(vertices[:,2]),np.min(vertices[:,2])])
            vertices[:,2] *= 0.5*ex.del_x[2]*ex.shape[2]/np.max(vertices[:,2])
            vertices[:,1] += (vertices[:,2] + 0.5*ex.del_x[2]*ex.shape[2])*ex.k0_vector[1]/ex.k0_vector[2]
            vertices[:,0] += (vertices[:,2] + 0.5*ex.del_x[2]*ex.shape[2])*ex.k0_vector[0]/ex.k0_vector[2]
            vertices *= scale
            # beam outside other end (exit)
            mesh  = drawing.add_isosurface_mesh(isosurface_mesh_from_matrix(input_beam, level = beam_level,
                                    step_size = beam_step_size, del_x = ex.del_x, meshcolor=[0,0,0,0], facecolor=[1,0.6,0,beam_opacity], abscolor=True))
            vertices = mesh.node_locations
            vertices[:,2] -= np.average([np.max(vertices[:,2]),np.min(vertices[:,2])])
            vertices[:,2] *= 0.5*ex.del_x[2]*ex.shape[2]/np.max(vertices[:,2])
            vertices[:,2] =  (vertices[:,2]*extend_beam[1]+(0.5+extend_beam[1]/2)*ex.del_x[2]*ex.shape[2]) # shift to back
            vertices[:,1] += (vertices[:,2] + 0.5*ex.del_x[2]*ex.shape[2])*ex.k0_vector[1]/ex.k0_vector[2]
            vertices[:,0] += (vertices[:,2] + 0.5*ex.del_x[2]*ex.shape[2])*ex.k0_vector[0]/ex.k0_vector[2]
            vertices *= scale
        else: # draw the beam as a plane line
            #z_vec = np.array([0,0,1])
            #beam_y = np.cross(ex.Q_vector,z_vec)
            #beam_y = beam_y/np.sqrt(np.sum(beam_y**2))
            node_locs = np.zeros((8,3))

            surface_angle = -np.arctan2(ex.k0_vector[1], -ex.k0_vector[0])
            mid_point = ex.params['BEAM']['mid'] # mid in relative units
            mid_point = [mid_point[0] * ex.shape[0]*ex.del_x[0], mid_point[1] * ex.shape[1]*ex.del_x[1]]

            thickness_offset = np.cross(np.array([-np.sin(surface_angle),np.cos(surface_angle),0]), ex.k0_vector)
            thickness_offset *= 0.2/np.sqrt(np.sum(thickness_offset**2))
            # front intesection of beam
            beam_loc_0 = np.array([0,1,0])*ex.params['BEAM']['transverse_width']
            node_locs[0,0] = np.cos(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]+beam_loc_0[0]+mid_point[0]) - np.sin(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]+beam_loc_0[1]+mid_point[1])
            node_locs[0,1] = np.cos(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]+beam_loc_0[1]+mid_point[1]) + np.sin(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]+beam_loc_0[0]+mid_point[0])
            node_locs[0,:] *= scale
            node_locs[1,0] = np.cos(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]-beam_loc_0[0] +mid_point[0]) - np.sin(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]-beam_loc_0[1] +mid_point[1])
            node_locs[1,1] = np.cos(surface_angle)*(-0.5*ex.shape[1]*ex.del_x[1]-beam_loc_0[1] +mid_point[1]) + np.sin(surface_angle)*(-0.5*ex.shape[0]*ex.del_x[0]-beam_loc_0[0] +mid_point[0])
            node_locs[1,:] *= scale

            node_locs[2,:] = node_locs[1,:]
            node_locs[3,:] = node_locs[0,:]
            # back intersection of beam
            node_locs[4,:] = node_locs[0,:] + box_shape[2]*ex.k0_vector/ex.k0_vector[2]
            node_locs[5,:] = node_locs[1,:] + box_shape[2]*ex.k0_vector/ex.k0_vector[2]

            node_locs[6,:] = node_locs[5,:]
            node_locs[7,:] = node_locs[4,:]
            # shift
            node_locs[:4,:] += -box_shape[2]*extend_beam[0]*ex.k0_vector[np.newaxis,:]/ex.k0_vector[2]
            node_locs[4:,:] += box_shape[2]*extend_beam[1]*ex.k0_vector[np.newaxis,:]/ex.k0_vector[2]
            node_locs[[0,1,4,5],:] += 0.5*thickness_offset
            node_locs[[2,3,6,7],:] -= 0.5*thickness_offset
            node_locs[:,:] += -0.5*box_shape[2]*np.array([0,0,1])

            box_node_collection = object_classes.NodeCollection(node_locs)
            drawing.nodes.append(box_node_collection)
            if 1:
                for seg in [[0,1,2,3],[4,5,6,7],[0,4,7,3],[7,3,2,6],[2,6,5,1],[5,1,0,4]]:
                    drawing.objects.append(object_classes.BoxFacet(box_node_collection, seg, facecolor = [1,0.6,0,beam_opacity], abscolor=False))

            #drawing.objects.append(object_classes.BoxFacet(box_node_collection, [0,1,2,3], facecolor = [1,0.6,0,beam_opacity*2], abscolor=False))
            #drawing.objects.append(object_classes.BoxFacet(box_node_collection, [0,1,2,3], facecolor = [1,0.6,0,beam_opacity*2], abscolor=False))


    # add internals
    print_time('# add internals')
    if not type(volume_data) == type(None):
        for mesh in volume_data:
            mesh.node_locations *= scale
            drawing.add_isosurface_mesh(mesh)
        '''vetices_list, faces_list = make_isosurface_and_get_vertice_collections(volume_data, level = volume_level, step_size = volume_step_size)
        for vertices, faces, i in zip(vetices_list, faces_list, range(len(faces_list))):
            vertices = vertices*ex.del_x-0.5*(ex.shape*ex.del_x) #scale to shape of box and center on zero
            drawing.add_mesh(vertices, faces, meshcolor=[0,0,0,0], facecolor=[0.8,0.2,0,0.1], abscolor=True)'''

    #lens
    if lens_scale>0:
        lens_radius = lens_scale*0.5*np.sqrt(2)*scale*np.max([ex.del_x[0]*ex.shape[0],ex.del_x[1]*ex.shape[1]])
        lens_pos = [0, 0, 0.5*scale*0.5*ex.del_x[2]*ex.shape[2]*(1+2*extend_imaging_system[1])]
        lens_pos[0] +=  (lens_pos[2]-0.5*box_shape[2])*ex.kh_vector[0]/ex.kh_vector[2]
        lens_pos[1] +=  (lens_pos[2]-0.5*box_shape[2])*ex.kh_vector[1]/ex.kh_vector[2]
        r_curvature = lens_radius*3

        delta_lens = ex.kh_vector/np.sqrt(np.sum(ex.kh_vector**2))
        delta_lens*= 1.1*2*(r_curvature-np.sqrt(r_curvature**2-lens_radius**2))
        num_lenses = 3
        for i in np.arange(num_lenses):
            i-= (num_lenses-1)/2-1
            drawing.add_lens(radius = lens_radius, num_links = 15, r_curvature = r_curvature,
                    facing = -ex.kh_vector, displacement = lens_pos+delta_lens*i)

    # add integrated image
    print_time('# add integrated image')
    if not type(intensity_imaging_system)==type(None):
        pcolormesh = object_classes.QuadrilateralColormesh(det_projection_nodes, [4,5,6,7],
                                    intensity_imaging_system, opacity_data = intensity_imaging_system_opacity, colormap = img_colormap)
        drawing.objects.append(pcolormesh)
    # add detector image
    '''if not type(detector_intensity)==type(None):
        pcolormesh = object_classes.QuadrilateralColormesh(det_detector_nodes, [4,5,6,7],
                                    detector_intensity, colormap = img_colormap)
        drawing.objects.append(pcolormesh)'''

    # add Q arrow
    if draw_scattering_arrows:
        Q_arrows = []
        Q_norm = ex.Q_vector/(np.sqrt(np.sum(ex.k0_vector**2))+np.sqrt(np.sum(ex.kh_vector**2)))
        ar = drawing.add_arrow([0,0,0], Q_norm*9,color=[0,0,0,1],nodes_per_circle=arrow_nodes_per_circle,head_fraction = 0.45/1, rod_radius = 0.1, hat_radius = 0.25)
        Q_arrows.append(ar)
        ar = drawing.add_text(Q_norm*9, 'Q')
        Q_arrows.append(ar)
        k0_norm = ex.k0_vector/np.sqrt(np.sum(ex.k0_vector**2))
        ar = drawing.add_arrow(-k0_norm*4.5, [0,0,0],color=[0.5,0,0,1],nodes_per_circle=arrow_nodes_per_circle,head_fraction = 0.45/3, rod_radius = 0.1/3, hat_radius = 0.25/3)
        k0_butt_node = ar.node_locations[0] # use this to rotate stage
        Q_arrows.append(ar)
        ar = drawing.add_text(-k0_norm*4.5+np.array([1,0,0]), 'k0')
        Q_arrows.append(ar)
        kh_norm = ex.kh_vector/np.sqrt(np.sum(ex.kh_vector**2))
        ar = drawing.add_arrow([0,0,0], kh_norm*4.5,color=[0.5,0,0,1],nodes_per_circle=arrow_nodes_per_circle,head_fraction = 0.45/3, rod_radius = 0.1/3, hat_radius = 0.25/3)
        Q_arrows.append(ar)
        ar = drawing.add_text(kh_norm*4.5, 'kh')
        Q_arrows.append(ar)


    ax_hkl = ['x','y','z']
    if 0:
        if not type(crystal_rotation_function) == type(None):
            '''
            we here calculat the vector of the sample suraces in the crystal basis
            i.e. the inverse of crystal_rotation_function
            (crystal_rotation_function gives a the rotation of the crystal vrt the sample geometry)
            '''
            rot_mat = np.stack((crystal_rotation_function(np.array([1.0,0,0])),
                                crystal_rotation_function(np.array([0,1.0,0])),
                                crystal_rotation_function(np.array([0,0,1.0])))).T
            rot_mat_inv =  np.linalg.inv(rot_mat)
            hkl_ax = vector_to_hkl(rot_mat_inv[:,0])
            ax_hkl[0] += f' = [{hkl_ax[0]},{hkl_ax[1]},{hkl_ax[2]}]'
            hkl_ax = vector_to_hkl(rot_mat_inv[:,1])
            ax_hkl[1] += f' = [{hkl_ax[0]},{hkl_ax[1]},{hkl_ax[2]}]'
            hkl_ax = vector_to_hkl(rot_mat_inv[:,2])
            ax_hkl[2] += f' = [{hkl_ax[0]},{hkl_ax[1]},{hkl_ax[2]}]'
    if draw_face_axes:
        ar = drawing.add_arrow([0,0,0], [3,0,0],color=[1,0.1,0,1],nodes_per_circle=arrow_nodes_per_circle, angle_offset=0)
        ar *= 30/8
        ar += np.array([1,0,0])*box_shape[0]*0.5
        ar = drawing.add_arrow([0,0,0], [0,3,0],color=[0,0.9,0.3,1],nodes_per_circle=arrow_nodes_per_circle)
        ar *= 30/8
        ar += np.array([0,1,0])*box_shape[1]*0.5
        ar = drawing.add_arrow([0,0,0], [0,0,3],color=[0,0.1,1,1],nodes_per_circle=arrow_nodes_per_circle)
        ar *= 30/8
        ar += np.array([box_shape[0]*0.25,box_shape[1]*0.25,0]) # shift this to off-center
        ar += np.array([0,0,1])*box_shape[2]*0.5
        ar = drawing.add_text([2.9,0,0], ax_hkl[0])#'[100]')
        ar *= 30/8
        ar += np.array([1,0,0])*box_shape[0]*0.5
        ar = drawing.add_text([-1.3,2.8,0.1], ax_hkl[1])#'[010]')
        ar *= 30/8
        ar += np.array([0,1,0])*box_shape[1]*0.5
        ar = drawing.add_text([0,0,2.9], ax_hkl[2])#'[001]')
        ar *= 30/8
        ar += np.array([0,0,1])*box_shape[2]*0.5
        ar += np.array([box_shape[0]*0.25,box_shape[1]*0.25,0]) # shift this to off-center

    # add xyz arrows
    if draw_axes:
        xyz_arrows = []
        print_time('# add arrows')
        ar = drawing.add_arrow([0,0,0], [0,0,3],color=[0,0.1,1,1],nodes_per_circle=arrow_nodes_per_circle)
        xyz_arrows.append(ar)
        ar = drawing.add_arrow([0,0,0], [0,3,0],color=[0,0.9,0.3,1],nodes_per_circle=arrow_nodes_per_circle)
        xyz_arrows.append(ar)
        ar = drawing.add_arrow([0,0,0], [3,0,0],color=[1,0.1,0,1],nodes_per_circle=arrow_nodes_per_circle, angle_offset=0)
        xyz_arrows.append(ar)

        ar = drawing.add_text([2.9,0,0], ax_hkl[0])#'[100]')
        xyz_arrows.append(ar)
        ar = drawing.add_text([-1.3,2.8,0.1], ax_hkl[1])#'[010]')
        xyz_arrows.append(ar)
        ar = drawing.add_text([0,0,2.9], ax_hkl[2])#'[001]')
        xyz_arrows.append(ar)

    # add curved arrows
    if draw_curved_arrows:
        curved_arrows = []
        # rocking curve
        curved_ar_links=21
        ar = drawing.add_curved_arrow([0,0,0], [5,0,0],[0,0,-1], angular_extent = np.pi/3, links = curved_ar_links,
                             head_fraction = 0.2, rod_radius = 0.05, hat_radius = 0.1, nodes_per_circle = arrow_nodes_per_circle, color=[1,0,1,1.0])
        curved_arrows.append(ar)
        ar = drawing.add_text(np.copy(ar.node_locations[-1]), 'Rocking')
        curved_arrows.append(ar)
    #drawing.add_isosurface_mesh(top_stage_mesh) # add top of stage

    # top part of stage
    if draw_stage:
        stage_dists = np.array([18.5,23,27,28,30,30.25,32.25, 32.5])-5
        stage_dists *= scale
        num_along_edge = 20
        stage_size = 20
        nodes, faces = add_stage_section(stage_dists[0], stage_dists[1], num_along_edge = num_along_edge, stage_size = stage_size)
        nodes[:4*num_along_edge,0] = -stage_dists[0]
        drawing.add_mesh(nodes, faces, facecolor=[0.4,0.3,0.2,1], abscolor=False)

    if not type(tilt_to_vector) == type(None): # tilt fig rocking
        drawing.rot_y(-np.arctan(tilt_to_vector[2]/tilt_to_vector[0]))
    # mid part of stage
    if draw_stage:
        nodes, faces = add_stage_section(stage_dists[1]+0.1, stage_dists[2], num_along_edge = num_along_edge, stage_size = stage_size)
        temp = np.copy(nodes[:,1])
        nodes[:,1] = nodes[:,2]
        nodes[:,2] = temp
        drawing.add_mesh(nodes, faces, facecolor=[0.5,0.4,0.3,1], abscolor=False)
    # azimuthal curve
    #drawing.add_isosurface_mesh(mid_stage_mesh) # add middle part of stage
    if draw_curved_arrows:
        ar = drawing.add_curved_arrow([0,0,0],[6,0,0],[0,1,0], angular_extent = np.pi/3, links = curved_ar_links,
                             head_fraction = 0.2, rod_radius = 0.05, hat_radius = 0.1, nodes_per_circle = arrow_nodes_per_circle, color=[0.7,0.7,0,1.0])
        curved_arrows.append(ar)
        ar = drawing.add_text(np.copy(ar.node_locations[-1]), 'Azimuthal')
        curved_arrows.append(ar)

    if not type(tilt_to_vector) == type(None): # tilt fig azimuthal
        drawing.rot_z(np.arctan2(tilt_to_vector[1],np.sqrt(tilt_to_vector[0]**2+tilt_to_vector[2]**2)))

    # bottom part of stage
    if draw_stage:
        nodes, faces = add_stage_section(stage_dists[2]+0.1, stage_dists[3], num_along_edge = num_along_edge, stage_size = stage_size)
        nodes[4*num_along_edge:,0] = -stage_dists[3]+0.1
        drawing.add_mesh(nodes, faces, facecolor=[0.6,0.5,0.4,1], abscolor=False)

        # stage rotator
        # top rotator
        ar = drawing.add_arrow([-stage_dists[3],0,0],
                               [-stage_dists[4],0,0],
                               head_fraction = 0, rod_radius = 0.8*stage_size/(stage_dists[4]-stage_dists[3]),
                               hat_radius = 0.1, nodes_per_circle = num_along_edge*4, color=[0.5,0.4,0.3,1], angle_offset=0)
        # handle dot
        ar = drawing.add_arrow([-0.5*(stage_dists[3]+stage_dists[4]),0.85*stage_size,0],
                               [-0.5*(stage_dists[3]+stage_dists[4]),0.8*stage_size,0],
                               head_fraction = 0, rod_radius = 17*0.3*(stage_dists[4]-stage_dists[3])/(0.55*stage_size),
                               hat_radius = 0.0001, nodes_per_circle = num_along_edge, color=[0.5,0.4,0.3,1], angle_offset=0)
        if not type(tilt_to_vector) == type(None):  # tilt fig rotation
            drawing.rot_x(-np.arcsin(k0_butt_node[1]/np.sqrt(k0_butt_node[1]**2+k0_butt_node[2]**2)))
        # bottom rotator
        ar = drawing.add_arrow([-stage_dists[5],0,0],
                               [-stage_dists[6],0,0],
                               head_fraction = 0, rod_radius = 0.8*stage_size/(stage_dists[6]-stage_dists[5]),
                               hat_radius = 0.1, nodes_per_circle = num_along_edge*4, color=[0.6,0.5,0.4,1], angle_offset=0)
        # handle dot
        ar = drawing.add_arrow([-0.5*(stage_dists[5]+stage_dists[6]),0.85*stage_size,0],
                               [-0.5*(stage_dists[5]+stage_dists[6]),0.8*stage_size,0],
                               head_fraction = 0, rod_radius = 17*0.3*(stage_dists[4]-stage_dists[3])/(0.55*stage_size),
                               hat_radius = 0.0001, nodes_per_circle = num_along_edge, color=[0.6,0.5,0.4,1], angle_offset=0)
        # table
        #table_nodes, _, __ = drawing.add_box([0,stage_size*10,stage_size*10], [-stage_dists[7],0,0], angle_vector=np.array([0,0,0]),
        #                                   meshcolor=[0,0,0,0], linewidth = 0, facecolor=[0.7,0.7,0.7,1], abscolor = False)
    # rotation curve
    if draw_curved_arrows:
        f=1.5/4
        ar = drawing.add_curved_arrow([-1,0,0],[-1,6/np.sqrt(2)/1.5,-6/np.sqrt(2)/1.5],[0,0,1], angular_extent = 4*np.pi/3, links = curved_ar_links*3,
                             head_fraction = 0.1, rod_radius = 0.05*f, hat_radius = 0.1*f, nodes_per_circle = arrow_nodes_per_circle, color=[0,0.7,0.7,1.0], angle_offset=np.pi)
        curved_arrows.append(ar)
        ar = drawing.add_text([-1,1,4], 'Rotation')
        curved_arrows.append(ar)
    # rotate perspective nicely
    print_time('# rotate perspective nicely')
    drawing.rot_y(-90*np.pi/180)
    drawing.rot_x(100*np.pi/180) # 95
    drawing.rot_y(-20*np.pi/180)
    # translate arrows
    if draw_scattering_arrows:
        for ar in Q_arrows:
            ar*=30/8
            ar += np.array([15,-10,10])*30/8
    if draw_axes:
        for ar in xyz_arrows:
            ar*=30/8
            ar += np.array([-5,-10,10])*30/8
    if draw_curved_arrows:
        for ar in curved_arrows:
            ar*=30/8
            ar += np.array([15,-10,10])*30/8
    # move crystal_structure
    if not type(atom_list) == type(None):
        for atom in crystal_structure:
            atom += crysta_structure_position
    # final rotation
    if not type(final_rotation_function) == type(None):
        final_rotation_function(drawing)


    # draw
    print_time('# draw')
    drawing.draw_structure(fig, ax, factor = factor, view_angle=view_angle, axis_off=True)
    print_time('#  rendered')
    # enable render
    if object_classes.backend == 'mayavi':
        fig.scene.disable_render = False
    if object_classes.backend == 'plotly':
        import plotly.graph_objects
        fig = plotly.graph_objects.Figure(data=ax['data'])

        camera = dict(
            up=dict(x=0, y=1, z=0),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=0, y=0, z=2),
            projection = dict(
            type="orthographic",
            ),
        )

        scene = dict(
                xaxis = dict(visible=False),
                yaxis = dict(visible=False),
                zaxis = dict(visible=False),
                annotations = ax['annotations'],
            )
        fig.update_layout( scene_camera = camera, scene = scene,
                margin=dict(t=0, r=0, l=0, b=0),
                height=800,
                width = 1100,
                )
    if not type(export_blender_filepath) == type(None):
        drawing.export_blender(export_blender_filepath)

    return drawing, fig, ax


def isosurface_mesh_from_matrix(threeD_matrix, level = 0, step_size = 10, del_x = np.array([1,1,1]), meshcolor=[0,0,0,0], facecolor=[1,0.6,0,0.7], abscolor=True):
    '''
    Makes an isosurface from threeD_matrix using skimage.measure.marching_cubes and generates a mesh centered on zero
    input:
        threeD_matrix: 3d numpy array of intensity
        level: level for the isosurface
        step_size: step size for the triangle mesh, see skimage.measure.marching_cubes
        del_x: numpy array of shape (3), step size in threeD_matrix
        meshcolor: numpy array of shape (4), color of the mesh lines
        facecolor: numpy array of shape (4), color of the mesh faces
        abscolor: bool, default True. If True, all facets are given same color, if False, coror is dependent angle of facet
    returns
        mesh: a object_classes.Mesh object
    '''
    # generate isosrface
    import skimage.measure
    verts, faces, normals, values = skimage.measure.marching_cubes(threeD_matrix, level, step_size=step_size)
    # shape to real shape
    verts *= del_x*(np.array(threeD_matrix.shape)/(np.array(threeD_matrix.shape)-1))
    verts -= 0.5*(threeD_matrix.shape*del_x)
    # generate mesh
    mesh = object_classes.Mesh(verts, faces, meshcolor=meshcolor, facecolor=facecolor, abscolor = abscolor)
    return mesh

def split_mesh_in_bodies(verts, faces):
    '''
    Splits meshes into independent bodies. This helps the projection for matplotlib.
    input:
        verts, numpy array of floats with shape (N,3) # locations in space
        faces, numpy array of ints with shape (M,3) # faces connecting three nodes
    returns:
        vetices_list
        faces_list
    '''
    # find independent bodies and label the vertices into bodies
    vert_groupings = np.arange(len(verts))
    changed = True
    while changed:
        changed = False
        for i, face in enumerate(faces):
            m = np.min([vert_groupings[face[0]],vert_groupings[face[1]],vert_groupings[face[2]]])
            if not all(vert_groupings[face]==m):
                vert_groupings[face]=m
                changed=True
    groups = list(set(vert_groupings))

    # make lists of faces by finding what group each face is in
    faces_list = []
    for i, group in enumerate(groups):
        faces_list.append([])
        for face in faces:
            if vert_groupings[face[0]]==group:
                faces_list[-1].append(face)

    # make lists of vertices and adjust the integers in faces_list to correspond to new numbers
    vetices_list = []
    for j, faces in enumerate(faces_list):
        faces_list[j] = np.array(faces)
        faces = faces_list[j]
        vertices_in_set = list(set(faces.ravel()))
        current_vertices = np.array(verts[vertices_in_set])
        vetices_list.append(current_vertices)
        tempmap = {}
        for i, vertice in enumerate(vertices_in_set):
            tempmap[vertice] = i
        for face in faces:
            face[0] = tempmap[face[0]]
            face[1] = tempmap[face[1]]
            face[2] = tempmap[face[2]]

    return vetices_list, faces_list



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

    def add_isosurface_mesh(self, mesh):
        '''
        adds a ready formed object_classes.Mesh
        input:
            mesh: object_classes.Mesh object
        returns:
            mesh: object_classes.Mesh object
        '''
        self.nodes.append(mesh)
        self.objects.append(mesh)
        return mesh

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

    def add_pcolormesh(self, node_locations,  intensity_imaging_system):
        '''
        adds a pcolormesh spanning node locations forming a quadrilateral
        input:
            node_locations: list of four node locations (np.array with shape (4,3))
            intensity_imaging_system: intensity to display on the pcolormesh
        returns:
            object_classes.NodeCollection object
            object_classes.QuadrilateralColormesh object
        '''
        node_collection = object_classes.NodeCollection(node_locations)
        self.nodes.append(box_node_collection)
        pcolormesh = object_classes.QuadrilateralColormesh(node_collection, [0,1,2,3], intensity_imaging_system)
        self.objects.append(pcolormes)
        return node_collection, pcolormesh

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


    def draw_structure(self, fig, ax, axis_off = True, xlim = None, ylim = None, factor = 1000, view_angle=8):
        '''
        Renders the scene
        input:
            fig: matplotlib figure
            ax: matplotlib axes
            xlim: optional default = None, the xlimits in the axis [xmin, xmax]
            ylim: optional default = None, the ylimits in the axis [ymin, ymax]
            axis_off: optional default = True, if True disables the axis fromain gthe scene
            axis_off: optional default = True, if True disables the axis fromain gthe scene
            axis_off: optional default = True, if True disables the axis fromain gthe scene
            factor: optional default = 1000, distance from object to camera
            view_angle: mayavi view angle, optional
        '''
        self.fig = fig
        self.ax = ax
        if object_classes.backend == 'matplotlib':
            # transform objects from 3d to 2d plane
            for object in self.objects:
                object.transform(factor=factor)
                if not hasattr(object,'layer'):
                    object.layer=0
            # sort objects based distance from center of object to camera
            self.objects=sorted(self.objects, key=lambda x: (x.loc[0]**2+x.loc[1]**2+(x.loc[2]-factor)**2), reverse=True)
            # sort objects based manually set labels, if any
            self.objects=sorted(self.objects, key=lambda x: x.layer, reverse=False)
            # add order as a parameter to the objects
            for i, object in enumerate(self.objects):
                object.z_order = i
        # draw objects
        for obj in self.objects:
            obj.draw(self.ax)
            print_time_render('# rendered '+str(obj))
        if object_classes.backend == 'matplotlib':
            # misc figure stuff
            if axis_off==True:
                self.ax.axis('off')
            self.ax.set_aspect('equal')
            if not type(xlim)==type(None):
                self.structure.ax.set_xlim(xlim)
            if not type(ylim)==type(None):
                self.structure.ax.set_xlim(ylim)
        elif object_classes.backend == 'mayavi':
            import mayavi.mlab
            mayavi.mlab.view(azimuth=0, elevation=0, distance=factor, focalpoint=[0,0,0], figure=fig)
            camera = fig.scene.camera
            camera.view_angle = view_angle


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
                    if 0: #circular
                        dist = i*branch_length
                        nodes2[-1] *= dist/np.sqrt(np.sum(node[0:3]**2))
                    elif 0: #square by cut
                        dist = i*branch_length
                        nodes2[-1] *= dist/np.sqrt(np.sum(node[0:3]**2))
                        if nodes2[-1][0] > sqare_lim: nodes2[-1][0]=sqare_lim
                        if nodes2[-1][0] < -sqare_lim: nodes2[-1][0]=-sqare_lim
                        if nodes2[-1][1] > sqare_lim: nodes2[-1][1]=sqare_lim
                        if nodes2[-1][1] < -sqare_lim: nodes2[-1][1]=-sqare_lim
                    else: #square move
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
            if 0:
                # if we convert it to a square by cut, we have to remove vertices that no longer fit
                for i in np.arange(len(facets))[::-1]:
                    facet = facets[i]
                    p0 = nodes[int(facet[0])]
                    p1 = nodes[int(facet[1])]
                    p2 = nodes[int(facet[2])]
                    normals = np.cross(p0 - p1, p2 - p1)
                    if sum(normals**2) < 10**-20:
                        facets.pop(i)
            facets = np.array(facets, dtype = int)
            # rotate from facing a to facing b
            facing /= np.sqrt(np.sum(facing**2))
            original_facing = np.array([0, 0, 1.0])
            if not all(facing == original_facing):
                #https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
                v = np.cross(original_facing, facing)
                c = np.dot(original_facing, facing)
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

        nodes += np.array(displacement)
        mesh = object_classes.Mesh(nodes, facets, meshcolor=[0,0,0,0], facecolor=[0.9,0.9,0.9,1], abscolor=False)
        self.nodes.append(mesh)
        self.objects.append(mesh)

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


def make_meshes_from_domains(domains, del_x, levels = None, colors = None, mod=1, step_size = 10, edge_zero = True):
    '''
    function that generates a list of dfxrm.three_d_draw_object_classes.Mesh objects from the domains object
    input:
        domains: 3d numpy array with a areas of constant vlaue signifying domains
        del_x: numpy array of set sizes for the domain array, i.e. ex.del_ex
        levels: optional default None. levels in domains to make meshes for. If None, a list is made by np.unique(domains)
        colors: optional, list of color [r,g,b,a] for different domanis. By default a list with 9 elements is used
        mod: optional float, default 1. Modulus to use for the domains array when assigning colors.
        step_size: int default 10, step size for the triangular mesh generation (skimage.measure.marching_cubes)
        edge_zero: bool default True. If true, will created mesh at edges in all cases
    '''
    def get_domains(domains, level, edge_zero):
        '''
        returns a 'boolean' copy of domains where everywhere (domains==level) will be set to 1, and all others are zero
        '''
        domains_to_show=np.copy(domains)
        domains_to_show[domains_to_show==level]=-1
        domains_to_show[domains_to_show!=-1]=0
        if edge_zero == True:
            # make outer edges zero:
            domains_to_show[0,:,:]=0
            domains_to_show[:,0,:]=0
            domains_to_show[:,:,0]=0
            # make sure a zero is within the step size at end
            domains_to_show[:,-step_size:,:]=0
            domains_to_show[-step_size:,:,:]=0
            domains_to_show[:,:,-step_size:]=0
        return -domains_to_show
    meshes=[]
    if type(colors)==type(None):
        colors = np.array([[1,0,0,0.95],
                           [0,1,0.3,0.95],
                           [0,0,1,0.95],
                           [0,1,1,0.95],
                           [1,0,1,0.95],
                           [1,1,0,0.95],
                           [0.25,0.25,0.25,0.95],
                           [0.5,0.5,.5,0.95],
                           [1,1,1,0.95],
                          ])
    if type(levels) == type(None):
        levels =  np.unique(domains)
    i = 0
    for dom_type in levels:
        try:
            if dom_type % mod == 0:
                i=0
            else:
                i+=1
            meshes.append(isosurface_mesh_from_matrix(get_domains(domains, dom_type, edge_zero), level = 0.9,
                                            step_size = step_size, del_x = del_x, meshcolor=[0,0,0,0],
                                            facecolor=colors[i%len(colors)], abscolor=False))
        except ValueError:
            print(f"ValueError in make_internals():get_domains()!  perhaps {dom_type} not in domains with step_size {step_size}.",
                  "You can safely ignore this message if it looks correct")
    return meshes

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
