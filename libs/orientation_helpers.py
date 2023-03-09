import numpy as np

def direction_from_string(xtl, sting):
    '''
    Takes as string '(h,k,l)' or '[u,v,w,]' or '' and converts to a real-space vector  
    input:
        xtl: dans_diffraction xtl object
        string: string of type '(h, k, l)' or [u, v, w]
    return:
        direction: numpy array of length 3
    '''
    s = sting.strip('[()]').split(',')
    arr = np.array([float(s[0]),float(s[1]),float(s[2])])
    if '(' in sting:
        direction = xtl.Cell.calculateQ(arr)[0]
    else:
        direction = xtl.Cell.calculateR(arr)[0]

    return direction

def Q_from_hkl_string(xtl, hkl_str):
    '''
    converts from a hkl value as string to a scattering vector
    input: 
        xtl: dans_diffraction xtl object
        hkl_str: string of type '(h, k, l)', i.e '(1,0,2)'
    return:
        hkl: hkl values, numpy array of lenght three 
        Q: real-space scattering vector, numpy array of lenght three 
    '''
    hkl_spli = hkl_str.strip('()').split(',')
    h = int(hkl_spli[0])
    k = int(hkl_spli[1])
    l = int(hkl_spli[2])
    hkl = np.array([h,k,l]) #'3, 3, 3'
    Q = xtl.Cell.calculateQ(hkl)[0]
    return hkl, Q
    
    
def up_dir_from_string(xtl, up_hkl_str, Q):
    
    '''
    Takes the sample 'up' as a string '(h,k,l)' or '[u,v,w,]' or '' and converts to a real-space vector  
    Note that orthogonality is enforced
    input: 
        xtl: dans_diffraction xtl object
        up_hkl_str: string of type '(h,k,l)' or '[u,v,w,]' or ''
        Q: numpy array of length 3, returned if up_hkl_str == ''
    return:
        up_dir: real-space vector, numpy array of lenght three 
    '''
    if not up_hkl_str == '':
        return direction_from_string(xtl, up_hkl_str)
    else:
        return Q

def front_dir_from_string(xtl, front_hkl_str, Q):
    
    '''
    Takes the sample 'front' as a string '(h,k,l)' or '[u,v,w,]' or '' and converts to a real-space vector  
    Note that orthogonality is enforced
    input: 
        xtl: dans_diffraction xtl object
        front_hkl_str: string of type '(h,k,l)' or '[u,v,w,]' or ''
        Q: numpy array of length 3, used if front_hkl_str = ''
    return:
        front_dir: real-space vector, numpy array of lenght three 
    '''
    if not front_hkl_str == '':
        front_dir = direction_from_string(xtl, front_hkl_str)
    else:
        front_dir = np.cross(Q, np.array([0,0,1]))
        if np.abs(np.sum(front_dir))<0.0001:
            front_dir = np.cross(Q,np.array([0,1,0]))
    return front_dir



########################## crystal rotation functions #####################
def rotate(x,z,phi):# 2D rotation function
    '''
    generic rotation function for generating the crystal rotation functions
    '''
    return x*np.cos(phi)-z*np.sin(phi), z*np.cos(phi)+x*np.sin(phi)

def rotate_x(loc,phi):#rotate_x: rotates around x axis (horisontal)
    '''
    in-place rotation around x axis for generating the crystal rotation function
    '''
    x,z=rotate(loc[1],loc[2],phi)
    loc[1] = x
    loc[2] = z
    
def rotate_y(loc,phi):# rotate_y: rotates around y axis (vertical)
    '''
    in-place rotation around y axis for generating the crystal rotation function
    '''
    x,z=rotate(loc[0],loc[2],phi)
    loc[0] = x
    loc[2] = z
    
def rotate_z(loc,phi):#rotate_z: rotates around z axis (out-of-plane)
    '''
    in-place rotation around z axis for generating the crystal rotation function
    '''
    x,z=rotate(loc[0],loc[1],phi)
    loc[0] = x
    loc[1] = z

class crystal_rotation:
    def __init__(self, z_rot, y_rot, x_rot):
        self.z_rot = z_rot
        self.y_rot = y_rot
        self.x_rot = x_rot
    
    def rotate(self, hkl):
        rotate_z(hkl, self.z_rot*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        rotate_y(hkl, self.y_rot*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        rotate_x(hkl, self.x_rot*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        return hkl
