import numpy as np
import copy
import plotly.graph_objects

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

brightness_change_from_angle = 0.5
light_direction = np.array([0.3, 1, 0.25])
light_direction = light_direction/np.sqrt(np.sum(light_direction**2))

def get_color(color, normal):
        '''
        Calculates the color based on the angle of the normal
        -> brightness with respect to the normal and the light_direction
        For simplicity, the object is as illuminated by both light_direction and -light_direction
        input:
            color: np.array of size 4 [r,g,b,a]
            normal: np.array of size 3 [x,y,z] or 2d np.array of shape (N,3)
        global values:
            brightness_change_from_angle: float from 0 to 1, fraction of change to apply
            light_direction: [x,y,z] the direction of illumination, must be normalized to length = 1
        output:
            color
        '''
        # ensure 2d input
        was_two_d = True
        if len(np.array(normal).shape) == 1: # 1D array:
            was_two_d = False
            normal = np.reshape(normal,(1,3))
        # prepare
        color=np.array(color, dtype = float)
        colors = np.zeros((normal.shape[0],4))
        # opacity # since we allow the user to move the figure in 3d we don't want to change the opacity
        '''thicknesses = 1/np.abs(normal[:,2])
        new_colors_3 = 1-(1-color[3])**thicknesses
        colors[:,3] = color[3] + (new_colors_3-color[3])*opacity_change_from_angle'''
        colors[:,3] = color[3]
        # brightness
        thickness_inv = np.abs(np.sum(normal*light_direction[np.newaxis,:], axis=1))
        f = (thickness_inv[thickness_inv>=0.5]-0.5)#*1.5
        colors[thickness_inv>=0.5,0:3] = (1-f[:,np.newaxis])*color[np.newaxis,0:3]+ f[:,np.newaxis]*np.ones((1,3))
        f = 1-(0.5-thickness_inv[thickness_inv<0.5])#*1.5
        colors[thickness_inv<0.5,0:3] = f[:,np.newaxis]*color[np.newaxis,0:3] #+ (1-f)*np.zeros((1,3))
        if was_two_d:
            return colors
        else:
            return colors[0]
class Text:
    def __init__(self, location, text, color=[0,0,0,1], scale =0.5):
        self.color=color
        self.loc = location
        self.text = text
        self.scale = scale

    def rotate_x(self,phi):
        '''
        rotate_x: rotates around x axis (horisontal)
        input: phi, angle in radians
        '''
        x,z=rotate(self.loc[1],self.loc[2],phi)
        self.loc[1] = x
        self.loc[2] = z
    def rotate_y(self,phi):
        '''
        rotate_y: rotates around y axis (vertical)
        input: phi, angle in radians
        '''
        x,z=rotate(self.loc[0],self.loc[2],phi)
        self.loc[0] = x
        self.loc[2] = z
    def rotate_z(self,phi):
        '''
        rotate_z: rotates around z axis (out-of-plane)
        input: phi, angle in radians
        '''
        x,z=rotate(self.loc[0],self.loc[1],phi)
        self.loc[0] = x
        self.loc[1] = z

    def draw(self,ax):
        '''
            Renders the text
        '''
        data = dict(
                showarrow=False,
                x=self.loc[0],
                y=self.loc[1],
                z=self.loc[2],
                text=self.text,
                xanchor="left",
                xshift=0,
                opacity=1)
        ax['annotations'].append(data)

    def __iadd__(self,vec):
        '''
        overload of the += operator to shift the location
        '''
        vec=np.array(vec)
        self.loc+=vec
        return self

    def __imul__(self,vec):
        '''
        overload of the *= operator to shift the location
        '''
        vec=np.array(vec)
        self.loc*=vec
        return self

    def to_dict(self):
        dic = {
            'type':'text',
            'location':self.loc,
            'text':self.text,
            'scale':self.scale,
        }
        return dic

class NodeCollection:
    '''
    Class that maintains node collections
        member variables:
            node_locations = locations of nodes in 3d
        member functions:
            rotate_x: rotates around x axis (horisontal)
            rotate_y: rotates around y axis (vertical)
            rotate_z: rotates around z axis (out-of-plane)
            +=: translation
            *=: magnify
    '''
    def __init__(self,node_locations):
        self.node_locations=np.array(node_locations)
    def rotate_x(self,phi):
        '''
        rotate_x: rotates around x axis (horisontal)
        input: phi, angle in radians
        '''
        for i, loc in enumerate(self.node_locations):
            x,z=rotate(loc[1],loc[2],phi)
            self.node_locations[i,1] = x
            self.node_locations[i,2] = z
    def rotate_y(self,phi):
        '''
        rotate_y: rotates around y axis (vertical)
        input: phi, angle in radians
        '''
        for i, loc in enumerate(self.node_locations):
            x,z=rotate(loc[0],loc[2],phi)
            self.node_locations[i,0] = x
            self.node_locations[i,2] = z
    def rotate_z(self,phi):
        '''
        rotate_z: rotates around z axis (out-of-plane)
        input: phi, angle in radians
        '''
        for i, loc in enumerate(self.node_locations):
            x,z=rotate(loc[0],loc[1],phi)
            self.node_locations[i,0] = x
            self.node_locations[i,1] = z
    def __iadd__(self,vec):
        '''
        overload of the += operator to shift the node locations
        '''
        vec=np.array(vec)
        self.node_locations+=vec[np.newaxis,:]
        return self

    def __imul__(self,vec):
        '''
        overload of the *= operator to shift the location
        '''
        vec=np.array(vec)
        self.node_locations*=vec
        return self


class Mesh(NodeCollection):
    '''
    A 3d traingle mesh object to render the mesh
    inherits from NodeCollection
    input:
        node_locations: numpy array of floats with shape (N,3) with coordinates in 3d for the nodes
        faces: numpy array of ints with shape (M,3) that make up the trangles
        meshcolor: float color of lines [r, g, b, a]
        facecolor: color of faces, [r, g, b, a]
    '''
    def __init__(self,node_locations, faces, meshcolor=[0,0,0,0], facecolor=[0,0.3,0.7,0.5], abscolor = False):
        self.meshcolor=meshcolor
        self.facecolor=facecolor
        self.faces = faces
        self.abscolor = abscolor
        #self.loc=np.average(self.node_locations, axis=0)
        #self.node_locations=np.array(node_locations)
        super().__init__(node_locations)
        return
    def get_colors(self):
        '''
        Get colors of each facet
        input:
            None
        return:
            colors: (N,4) matix corresponding to color of each facet
        '''
        colors = np.ones((self.faces.shape[0],4))
        p0 = self.node_locations[self.faces[:,0]]
        p1 = self.node_locations[self.faces[:,1]]
        p2 = self.node_locations[self.faces[:,2]]
        normals = np.cross(p0 - p1, p2 - p1)
        normals[:,:] /= np.sqrt(np.sum(normals**2,axis=1))[:,np.newaxis]
        colors = get_color(self.facecolor, normals)
        return colors

    def draw(self,ax):
        '''
            adds this object to the 'data' entry in the plotly dictionary ax
        '''
        if not self.abscolor:
            colors = self.get_colors()
            data = plotly.graph_objects.Mesh3d(
                            x = self.node_locations[:,0],
                            y = self.node_locations[:,1],
                            z = self.node_locations[:,2],
                            facecolor = colors,
                            alphahull=0,
                            flatshading = True,
                            i = self.faces[:,0],
                            j = self.faces[:,1],
                            k = self.faces[:,2],
                          hoverinfo='skip')
            ax['data'].append(data)
        else:
            self.facecolor = np.array(self.facecolor)
            facecolor = np.array(self.facecolor[:3]*256, dtype = int)
            data = plotly.graph_objects.Mesh3d(
                            x = self.node_locations[:,0],
                            y = self.node_locations[:,1],
                            z = self.node_locations[:,2],
                            color = 'rgb'+str(tuple(facecolor)),
                            opacity = self.facecolor[3],
                            alphahull=0,
                            flatshading = True,
                            i = self.faces[:,0],
                            j = self.faces[:,1],
                            k = self.faces[:,2],
                          hoverinfo='skip')
            ax['data'].append(data)

    def to_dict(self):
        '''
        converts the object to a dictionary so that it can be exported to blender as part of a picle
        return:
            dictionary
        '''
        dic = {
            'type':'Mesh',
            'nodes':self.node_locations,
            'faces':self.faces,
            'meshcolor':self.meshcolor,
            'facecolor':self.facecolor,
            'abscolor':self.abscolor,
            'color':self.get_colors(),
        }
        return dic


        #quatplot(x,y, np.asarray(faces), ax=ax, )

class BoxFacet:
    '''
    A facet, used for example to make a box:
    input:
        node_collection: NodeCollection object, must cotain the nodes referenced in corner_nodes
        corner_nodes: nodes that make up the facet
        facecolor: color of faces, [r, g, b, a]
        abscolor: bool
            True -> all faces will be facecolor
            False -> opacity and color of face modified according to face normal with respect to camera
    '''
    def __init__(self, node_collection, corner_nodes, facecolor=[0,0.3,0.7,0.5], abscolor = False):
        self.node_collection = node_collection
        self.corner_nodes = corner_nodes
        self.facecolor = facecolor
        self.abscolor = abscolor

    def get_color(self):
        '''
        Get colors of the facet
        input:
            None
        return:
            colors: len (4) array corresponding to color + opacity
        '''
        normal = np.cross(   self.node_collection.node_locations[self.corner_nodes[0]] - self.node_collection.node_locations[self.corner_nodes[1]]
                            ,self.node_collection.node_locations[self.corner_nodes[2]] - self.node_collection.node_locations[self.corner_nodes[1]])
        normal = normal/np.sqrt(np.sum(normal**2))
        color = get_color(self.facecolor, normal)
        return color
    def draw(self,ax):
        '''
            adds this object to the 'data' entry in the plotly dictionary ax
        '''
        self.facecolor = np.array(self.facecolor)
        if not self.abscolor:
            colors = self.get_color()
        else:
            colors = self.facecolor
        #colors[:3] = np.array(colors[:3]*256,dtype = int)
        #colors *= 0.5
        #print(colors)
        n = self.node_collection.node_locations[self.corner_nodes]
        data = plotly.graph_objects.Mesh3d(
                        x = n[:,0],
                        y = n[:,1],
                        z = n[:,2],
                        facecolor = [colors,colors],
                        alphahull = 0,
                        flatshading = True,
                        i = [0,0],
                        j = [1,3],
                        k = [2,2],
                        hoverinfo='skip',
                        )
        ax['data'].append(data)


    def to_dict(self):
        '''
        converts the object to a dictionary so that it can be exported to blender as part of a picle
        return:
            dictionary
        '''
        dic = {
            'type':'BoxFacet',
            'span':self.node_collection.node_locations[self.corner_nodes],
            'facecolor':self.facecolor,
            'abscolor':self.abscolor,
            'color':self.get_color(),
        }
        return dic

class BoxLine:
    '''
    A line, used for example to make a box:
    input:
        node_collection: NodeCollection object, must cotain the nodes referenced in corner_nodes
        corner_nodes: nodes that make up the line
        meshcolor: float color of lines [r, g, b, a]
        linewidth: float linewidth
    '''
    def __init__(self, node_collection, corner_nodes, meshcolor=[0,0,0,0], linewidth = 0.1):
        self.node_collection = node_collection
        self.corner_nodes = corner_nodes
        self.meshcolor = meshcolor
        self.linewidth = linewidth
    def draw(self,ax):
        '''
            adds this object to the 'data' entry in the plotly dictionary ax
        '''
        locs = self.node_collection.node_locations[self.corner_nodes]
        data = plotly.graph_objects.Scatter3d(
                        x = locs[:,0],
                        y = locs[:,1],
                        z = locs[:,2],
                        mode = "lines",
                        line = dict(
                            color = 'rgb'+str(tuple(np.array(self.meshcolor)*256))
                        ),
                        #linsecolor = 'rgb'+str((0,0,0)),
                        showlegend = False,
                      hoverinfo='skip')
        ax['data'].append(data)

    def to_dict(self):
        '''
        converts the object to a dictionary so that it can be exported to blender as part of a picle
        return:
            dictionary
        '''
        dic = {
            'type':'BoxLine',
            'span':self.node_collection.node_locations[self.corner_nodes],
            'meshcolor':self.meshcolor,
            'linewidth':self.linewidth,
        }
        return dic



class Atom:
    '''
    Atom, used to draw the crystal structure:
    input:
        loc: position
        typ: type of atom, i.e. 'Ba'
        color: None of color[r, g b]
        radius: None of float
    '''
    def __init__(self, loc, typ, color = None, radius = None):
        self.loc = np.array(loc)
        self.typ = typ
        if type(color) == type(None):
            self.color = atom_color_dict[typ]
        else:
            self.color=color
        if type(radius) == type(None):
            self.radius = atomic_radius_dict[typ]
        else:
            self.radius = radius
        return
    def __iadd__(self,vec):
        '''
        overload of the += operator to shift the location
        '''
        vec = np.array(vec)
        self.loc += vec
        return self

    def __imul__(self,vec):
        '''
        overload of the *= operator to shift the location and increase the radius
        '''
        self.loc *= vec
        self.radius *= vec
        return self


    def draw(self,ax):
        '''
            adds this object as an isochahedron to the 'data' entry in the plotly dictionary ax
        '''
        vertices = isochahedron_vertices*self.radius+self.loc
        # get colors (same as for meshes)
        colors = np.ones((isochahedron_faces.shape[0],4))
        p0 = vertices[isochahedron_faces[:,0]]
        p1 = vertices[isochahedron_faces[:,1]]
        p2 = vertices[isochahedron_faces[:,2]]
        normals = np.cross(p0 - p1, p2 - p1)
        normals[:,:] /= np.sqrt(np.sum(normals**2,axis=1))[:,np.newaxis]
        colors = get_color(self.color, normals)
        data = plotly.graph_objects.Mesh3d(
                        x = vertices[:,0],
                        y = vertices[:,1],
                        z = vertices[:,2],
                        facecolor = colors,
                        alphahull=0,
                        flatshading = True,
                        i = isochahedron_faces[:,0],
                        j = isochahedron_faces[:,1],
                        k = isochahedron_faces[:,2],
                      hoverinfo='skip')
        ax['data'].append(data)


    def rotate_y(self,phi):
        '''
        rotate_y: rotates around y axis (vertical)
        input: phi, angle in radians
        '''
        x,z=rotate(self.loc[0],self.loc[2],phi)
        self.loc[0]=x
        self.loc[2]=z
    def rotate_x(self,phi):
        '''
        rotate_x: rotates around x axis (horizontal)
        input: phi, angle in radians
        '''
        x,z=rotate(self.loc[1],self.loc[2],phi)
        self.loc[1]=x
        self.loc[2]=z
    def rotate_z(self,phi):
        '''
        rotate_z: rotates around z axis (out-of-plane)
        input: phi, angle in radians
        '''
        x,z=rotate(self.loc[0],self.loc[1],phi)
        self.loc[0]=x
        self.loc[1]=z

    def to_dict(self):
        '''
        converts the object to a dictionary so that it can be exported to blender as part of a picle
        return:
            dictionary
        '''
        vertices = isochahedron_vertices*self.radius+self.loc
        # get colors colors (same as for meshes)
        colors = np.ones((isochahedron_faces.shape[0],4))
        p0 = vertices[isochahedron_faces[:,0]]
        p1 = vertices[isochahedron_faces[:,1]]
        p2 = vertices[isochahedron_faces[:,2]]
        normals = np.cross(p0 - p1, p2 - p1)
        normals[:,:] /= np.sqrt(np.sum(normals**2,axis=1))[:,np.newaxis]
        colors = get_color(self.color, normals)

        dic = {
            'type':'Mesh',
            'nodes':vertices,
            'faces':isochahedron_faces,
            'meshcolor':[0,0,0,0],
            'facecolor':self.color,
            'abscolor':False,
            'color':self.color,
        }
        return dic


radius_and_color_list=[
[1, 'H', 0.46, 1.20, 0.200, 1.00000, 0.80000, 0.80000],
[1, 'D', 0.46, 1.20, 0.200, 0.80000, 0.80000, 1.00000],
[2, 'He', 1.22, 1.40, 1.220, 0.98907, 0.91312, 0.81091],
[3, 'Li', 1.57, 1.40, 0.590, 0.52731, 0.87953, 0.45670],# 0.52731, 0.87953, 0.45670
[4, 'Be', 1.12, 1.40, 0.270, 0.37147, 0.84590, 0.48292],
[5, 'B', 0.81, 1.40, 0.110, 0.12490, 0.63612, 0.05948],
[6, 'C', 0.77, 1.70, 0.150, 0.50430, 0.28659, 0.16236],
[7, 'N', 0.74, 1.55, 1.460, 0.69139, 0.72934, 0.90280],
[8, 'O', 0.74, 1.52, 1.400, 0.99997, 0.01328, 0.00000],
[9, 'F', 0.72, 1.47, 1.330, 0.69139, 0.72934, 0.90280],
[10, 'Ne', 1.60, 1.54, 1.600, 0.99954, 0.21788, 0.71035],
[11, 'Na', 1.91, 1.54, 1.020, 0.97955, 0.86618, 0.23787],
[12, 'Mg', 1.60, 1.54, 0.720, 0.98773, 0.48452, 0.08470],
[13, 'Al', 1.43, 1.54, 0.390, 0.50718, 0.70056, 0.84062],
[14, 'Si', 1.18, 2.10, 0.260, 0.10596, 0.23226, 0.98096],
[15, 'P', 1.10, 1.80, 0.170, 0.75557, 0.61256, 0.76425],
[16, 'S', 1.04, 1.80, 1.840, 1.00000, 0.98071, 0.00000],
[17, 'Cl', 0.99, 1.75, 1.810, 0.19583, 0.98828, 0.01167],
[18, 'Ar', 1.92, 1.88, 1.920, 0.81349, 0.99731, 0.77075],
[19, 'K', 2.35, 1.88, 1.510, 0.63255, 0.13281, 0.96858],
[20, 'Ca', 1.97, 1.88, 1.120, 0.35642, 0.58863, 0.74498],
[21, 'Sc', 1.64, 1.88, 0.745, 0.71209, 0.38930, 0.67279],
[22, 'Ti', 1.47, 1.88, 0.605, 0.47237, 0.79393, 1.00000],
[23, 'V', 1.35, 1.88, 0.580, 0.90000, 0.10000, 0.00000],
[24, 'Cr', 1.29, 1.88, 0.615, 0.00000, 0.00000, 0.62000],
[25, 'Mn', 1.37, 1.88, 0.830, 0.66148, 0.03412, 0.62036],
[26, 'Fe', 1.26, 1.88, 0.780, 0.71051, 0.44662, 0.00136],
[27, 'Co', 1.25, 1.88, 0.745, 0.00000, 0.00000, 0.68666],
[28, 'Ni', 1.25, 1.88, 0.690, 0.72032, 0.73631, 0.74339],
[29, 'Cu', 1.28, 1.88, 0.730, 0.13390, 0.28022, 0.86606],
[30, 'Zn', 1.37, 1.88, 0.740, 0.56123, 0.56445, 0.50799],
[31, 'Ga', 1.53, 1.88, 0.620, 0.62292, 0.89293, 0.45486],
[32, 'Ge', 1.22, 1.88, 0.530, 0.49557, 0.43499, 0.65193],
[33, 'As', 1.21, 1.85, 0.335, 0.45814, 0.81694, 0.34249],
[34, 'Se', 1.04, 1.90, 1.980, 0.60420, 0.93874, 0.06122],
[35, 'Br', 1.14, 1.85, 1.960, 0.49645, 0.19333, 0.01076],
[36, 'Kr', 1.98, 2.02, 1.980, 0.98102, 0.75805, 0.95413],
[37, 'Rb', 2.50, 2.02, 1.610, 1.00000, 0.00000, 0.60000],
[38, 'Sr', 2.15, 2.02, 1.260, 0.00000, 1.00000, 0.15259],
[39, 'Y', 1.82, 2.02, 1.019, 0.40259, 0.59739, 0.55813],
[40, 'Zr', 1.60, 2.02, 0.720, 0.00000, 1.00000, 0.00000],
[41, 'Nb', 1.47, 2.02, 0.640, 0.29992, 0.70007, 1],# 0.29992, 0.70007, 0.46459
[42, 'Mo', 1.40, 2.02, 0.590, 0.70584, 0.52602, 0.68925],
[43, 'Tc', 1.35, 2.02, 0.560, 0.80574, 0.68699, 0.79478],
[44, 'Ru', 1.34, 2.02, 0.620, 0.81184, 0.72113, 0.68089],
[45, 'Rh', 1.34, 2.02, 0.665, 0.80748, 0.82205, 0.67068],
[46, 'Pd', 1.37, 2.02, 0.860, 0.75978, 0.76818, 0.72454],
[47, 'Ag', 1.44, 2.02, 1.150, 0.72032, 0.73631, 0.74339],
[48, 'Cd', 1.52, 2.02, 0.950, 0.95145, 0.12102, 0.86354],
[49, 'In', 1.67, 2.02, 0.800, 0.84378, 0.50401, 0.73483],
[50, 'Sn', 1.58, 2.02, 0.690, 0.60764, 0.56052, 0.72926],
[51, 'Sb', 1.41, 2.00, 0.760, 0.84627, 0.51498, 0.31315],
[52, 'Te', 1.37, 2.06, 2.210, 0.67958, 0.63586, 0.32038],
[53, 'I', 1.33, 1.98, 2.200, 0.55914, 0.12200, 0.54453],
[54, 'Xe', 2.18, 2.16, 0.480, 0.60662, 0.63218, 0.97305],
[55, 'Cs', 2.72, 2.16, 1.740, 0.05872, 0.99922, 0.72578],
[56, 'Ba', 2.24, 2.16, 1.420, 0.11835, 0.93959, 0.17565],
[57, 'La', 1.88, 2.16, 1.160, 0.35340, 0.77057, 0.28737],
[58, 'Ce', 1.82, 2.16, 0.970, 0.82055, 0.99071, 0.02374],
[59, 'Pr', 1.82, 2.16, 1.126, 0.99130, 0.88559, 0.02315],
[60, 'Nd', 1.82, 2.16, 1.109, 0.98701, 0.55560, 0.02744],
[61, 'Pm', 1.81, 2.16, 1.093, 0.00000, 0.00000, 0.96000],
[62, 'Sm', 1.81, 2.16, 1.270, 0.99042, 0.02403, 0.49195],
[63, 'Eu', 2.06, 2.16, 1.066, 0.98367, 0.03078, 0.83615],
[64, 'Gd', 1.79, 2.16, 1.053, 0.75325, 0.01445, 1.00000],
[65, 'Tb', 1.77, 2.16, 1.040, 0.44315, 0.01663, 0.99782],
[66, 'Dy', 1.77, 2.16, 1.027, 0.19390, 0.02374, 0.99071],
[67, 'Ho', 1.76, 2.16, 1.015, 0.02837, 0.25876, 0.98608],
[68, 'Er', 1.75, 2.16, 1.004, 0.28688, 0.45071, 0.23043],
[69, 'Tm', 1.00, 2.16, 0.994, 0.00000, 0.00000, 0.88000],
[70, 'Yb', 1.94, 2.16, 0.985, 0.15323, 0.99165, 0.95836],
[71, 'Lu', 1.72, 2.16, 0.977, 0.15097, 0.99391, 0.71032],
[72, 'Hf', 1.59, 2.16, 0.710, 0.70704, 0.70552, 0.35090],
[73, 'Ta', 1.47, 2.16, 0.640, 0.71952, 0.60694, 0.33841],
[74, 'W', 1.41, 2.16, 0.600, 0.55616, 0.54257, 0.50178],
[75, 'Re', 1.37, 2.16, 0.530, 0.70294, 0.69401, 0.55789],
[76, 'Os', 1.35, 2.16, 0.630, 0.78703, 0.69512, 0.47379],
[77, 'Ir', 1.36, 2.16, 0.625, 0.78975, 0.81033, 0.45049],
[78, 'Pt', 1.39, 2.16, 0.625, 0.79997, 0.77511, 0.75068],
[79, 'Au', 1.44, 2.16, 1.370, 0.99628, 0.70149, 0.22106],
[80, 'Hg', 1.55, 2.16, 1.020, 0.82940, 0.72125, 0.79823],
[81, 'Tl', 1.71, 2.16, 0.885, 0.58798, 0.53854, 0.42649],
[82, 'Pb', 1.75, 2.16, 1.190, 0.32386, 0.32592, 0.35729],
[83, 'Bi', 1.82, 2.16, 1.030, 0.82428, 0.18732, 0.97211],
[84, 'Po', 1.77, 2.16, 0.940, 0.00000, 0.00000, 1.00000],
[85, 'At', 0.62, 2.16, 0.620, 0.00000, 0.00000, 1.00000],
[86, 'Rn', 0.80, 2.16, 0.800, 1.00000, 1.00000, 0.00000],
[87, 'Fr', 1.00, 2.16, 1.800, 0.00000, 0.00000, 0.00000],
[88, 'Ra', 2.35, 2.16, 1.480, 0.42959, 0.66659, 0.34786],
[89, 'Ac', 2.03, 2.16, 1.120, 0.39344, 0.62101, 0.45034],
[90, 'Th', 1.80, 2.16, 1.050, 0.14893, 0.99596, 0.47106],
[91, 'Pa', 1.63, 2.16, 0.780, 0.16101, 0.98387, 0.20855],
[92, 'U', 1.56, 2.16, 0.730, 0.47774, 0.63362, 0.66714],
[93, 'Np', 1.56, 2.16, 0.750, 0.30000, 0.30000, 0.30000],
[94, 'Pu', 1.64, 2.16, 0.860, 0.30000, 0.30000, 0.30000],
[95, 'Am', 1.73, 2.16, 0.975, 0.30000, 0.30000, 0.30000],
[96, 'XX', 0.80, 1.00, 0.800, 0.30000, 0.30000, 0.30000]]
atomic_radius_dict={}
vdw_radius_dict={}
ionic_radius_dict={}
atom_color_dict={}
for e in radius_and_color_list:
    typ=e[1]
    atomic_radius_dict[typ]=e[2]/3
    vdw_radius_dict[typ]=e[3]/3
    ionic_radius_dict[typ]=e[4]/3
    atom_color_dict[typ]=np.array([e[5],e[6],e[7],1])

# isochahedron
phi = (1+np.sqrt(5))/2
isochahedron_vertices = np.array([
    [0,1,phi],  # 0
    [0,1,-phi], # 1
    [0,-1,-phi],# 2
    [0,-1,phi], # 3
    [phi,0,1],  # 4
    [-phi,0,1], # 5
    [-phi,0,-1],# 6
    [phi,0,-1], # 7
    [1,phi,0],  # 8
    [1,-phi,0], # 9
    [-1,-phi,0],# 10
    [-1,phi,0], # 11
    ])
isochahedron_vertices = isochahedron_vertices/np.sqrt(np.sum(isochahedron_vertices[0]**2))
for i, loc in enumerate(isochahedron_vertices):
    x,z=rotate(loc[0],loc[2],np.pi/5)
    isochahedron_vertices[i,0] = x
    isochahedron_vertices[i,2] = z
for i, loc in enumerate(isochahedron_vertices):
    y,z=rotate(loc[1],loc[2],0.03)
    isochahedron_vertices[i,1] = y
    isochahedron_vertices[i,2] = z

isochahedron_faces=np.array([
    [0,4,8], # all pos
    [1,7,8], # pos, pos, neg
    [3,4,9], # pos, neg, pos
    [0,5,11], # neg, pos, pos
    [2,7,9], # pos, neg, neg#
    [1,6,11], # neg, pos, neg#
    [3,5,10], # neg, neg, pos#
    [2,6,10], # neg, neg, neg
    [0,3,4],
    [0,3,5],
    [1,2,6],
    [1,2,7],#
    [4,7,8],
    [4,7,9],
    [5,6,10],
    [5,6,11],#
    [8,11,0],
    [8,11,1],
    [9,10,2],
    [9,10,3],#
])
