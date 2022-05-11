import numpy as np
import matplotlib.pyplot as plt
import Dans_Diffraction
import xraydb
def att_length(xtl,energies_kev):
    '''
    calculates attenuation length
    input:
        xtl, Dans_Diffraction.Crystal() object
        energies_kev, 1d array of floats or float, energies in keV
    return:
        1d numpy array of floats or float, attenuation length(s) in µm
    '''
    atom_type = xtl.Structure.type
    occ = xtl.Structure.occupancy
    natoms = np.sum(occ)
    vol = xtl.Properties.volume()*10**6/10**30 #Å^3
    weight = Dans_Diffraction.fc.atom_properties(atom_type, 'Weight')/(6.02214076*10**23)
    #st.write(weight)
    mu =  occ[0]*weight[0]*xraydb.mu_elam(atom_type[0], energies_kev*1000)
    #mu +=  occ[0]*weight[0]*xraydb.incoherent_cross_section_elam(atom_type[0], energies_kev*1000)
    # Values returned are in units of cm^2/gr # here mutiplied by weight
    for i in range(1,len(occ)):
        mu +=  occ[i]*weight[i]*xraydb.mu_elam(atom_type[i], energies_kev*1000)#photoabsorption_crosssection(elements, energy_kev)
        #mu +=  occ[0]*weight[0]*xraydb.incoherent_cross_section_elam(atom_type[0], energies_kev*1000)
    #st.write(vol)
    absorption = vol/mu
    return absorption*10**4 # cm to µm

def get_att_plot(xtl, energy_kev, sample_thickness):
    '''
    makes a plot of the attenuation length
        input:
            xtl, Dans_Diffraction.Crystal() object
            energy_kev, float, energy in keV
            sample_thickness, float, thickness of the sample in µm
        return:
            matplotlib figure
    '''
    fig, axes = plt.subplots(1,2,figsize=(12,4), dpi=100)

    energies = np.arange(1,30,0.1)
    attenuation_um2 = att_length(xtl, energies)

    axes[0].loglog(energies, attenuation_um2, color = [0.8,0.2,0])
    axes[0].set_xlabel('energy [keV]')
    axes[0].set_ylabel(r'Attenuation length [µm]')
    axes[0].set_xticks([1,10])
    axes[0].set_xticklabels([1,10])
    y_lim = axes[0].get_ylim()
    axes[0].plot([energy_kev, energy_kev], y_lim, color = [0,0.3,1,0.5])
    axes[0].set_ylim(y_lim)
    axes[0].text(energy_kev,y_lim[0]*2, f'E = {energy_kev:.2f} keV', rotation = 'vertical', ha = 'right')

    att = att_length(xtl, xtl.Scatter._energy_kev)
    axes[0].text(energy_kev,att, f'{att:.0f} µm', ha = 'right')

    thickness = np.arange(50,1000,1)

    xray_transmission = np.exp(-thickness/att)
    axes[1].loglog(thickness, xray_transmission, color = [0.8,0.2,0])
    axes[1].set_xlabel('Thickness [µm]')
    axes[1].set_ylabel(r'Transmission [-]')
    axes[1].set_title(f'Transmission @ {np.round(xtl.Scatter._energy_kev):.0f} keV')

    y_lim = axes[1].get_ylim()
    if y_lim[1]<10*y_lim[0]:
        y_lim=[0.1,1.1]
    axes[1].plot([sample_thickness, sample_thickness], y_lim, color = [0,0.3,1,0.5])
    axes[1].set_ylim(y_lim)
    axes[1].text(sample_thickness,y_lim[0]*2, f'Thickness = {sample_thickness:.2f} µm', rotation = 'vertical', ha = 'right')
    xray_transmission_0 = np.exp(-sample_thickness/att)
    axes[1].text(sample_thickness,xray_transmission_0, f'{xray_transmission_0*100:.2f}%')
    return fig

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
