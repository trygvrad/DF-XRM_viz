import streamlit as st
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import plotly
import copy

def plot_scatter(xtl, ax, energy_kev):
    uvw,atom_types,label,occ,uiso,mxmymz = xtl.Structure.get()
    a_scatter = xtl.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]
    pos = np.array([0.0,0.0])
    import three_d_draw_object_classes
    for scatter, atom, label in zip(a_scatter,atom_types,label):
        color = three_d_draw_object_classes.atom_color_dict[atom]
        #color = [float(color[0]),float(color[1]),float(color[2])]
        #ax.plot([pos[0],pos[0]+scatter.real],[pos[1],pos[1]+scatter.imag], color = [0,0,0], linewidth = 2.5)
        ax.plot([pos[0],pos[0]+scatter.real],[pos[1],pos[1]+scatter.imag],'--', color = color)
        tx = ax.text(pos[0]+0.5*scatter.real,pos[1]+0.5*scatter.imag, atom, color= color)
        tx.set_bbox(dict(facecolor=[1,1,1], alpha=0.8, edgecolor='None', pad=-0.8))
        pos[0] += scatter.real
        pos[1] += scatter.imag
    FS = np.sum(a_scatter)
    ax.plot([0,FS.real],[0,FS.imag], color = [0,0,0,1])
    # make arrowhead
    phi = np.arctan2(FS.imag,FS.real)
    L = np.sqrt(FS.imag**2+FS.real**2)*0.1
    FS_a0 = 1j*np.sin(phi+0.3)*L+np.cos(phi+0.3)*L
    FS_a1 = 1j*np.sin(phi-0.3)*L+np.cos(phi-0.3)*L
    ax.plot([FS.real-FS_a1.real, FS.real, FS.real-FS_a0.real],
            [FS.imag-FS_a1.imag, FS.imag, FS.imag-FS_a0.imag],
            color = [0,0,0,1])
    tx = ax.text(0.5*FS.real,0.5*FS.imag,r'$F$', color= [0,0,0])
    tx.set_bbox(dict(facecolor=[1,1,1], alpha=0.8, edgecolor='None', pad=-1.0))

    tx = ax.text(0,5*L,r'$|F|^2 = $'+f'{(FS.real**2+FS.imag**2):.6e}'+'\n'+r'$\phi_F = $'+f'{np.arctan2(FS.imag,FS.real)*180/np.pi:.6f}',
        ha = 'center',color= [0,0,0])
    tx.set_bbox(dict(facecolor=[1,1,1], alpha=0.8, edgecolor='None', pad=-1.0))

    ax.set_aspect(1)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    mx = np.max([-xlim[0], xlim[1], -ylim[0],ylim[1]])
    ax.plot([0,0], [-mx,mx], color = [0,0,0,1], lw = 0.1)
    ax.plot([-mx,mx],[0,0], color = [0,0,0,1], lw = 0.1)
    ax.set_xlim([-mx,mx])
    ax.set_ylim([-mx,mx])
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_xticks())
    ax.set_xlabel('$Re(F)$')
    ax.set_ylabel('$Im(F)$')
    ticks = ax.get_xticks()
    #ticks = ticks[(len(ticks)+1)//2:]
    ticks = np.arange(0,ticks[-1],1*(ticks[-1]-ticks[-2]))
    for tick in ticks:
        circle = plt.Circle((0, 0), tick, color=[0,0,0,0], ec= [0,0,0,1], lw = 0.1)
        ax.add_patch(circle)
    ax.set_xlim([-mx,mx])
    ax.set_ylim([-mx,mx])

st.markdown(
        f"""
<style>
    .reportview-container .main .block-container{{
        max-width: {1000}px;
        padding-top: {10}rem;
        padding-right: {2}rem;
        padding-left: {2}rem;
        padding-bottom: {10}rem;
    }}
    .reportview-container .main {{
        color: {"#000000"};
        background-color: {"#EAEAEA"};
    }}
</style>
""",
        unsafe_allow_html=True,
    )

st.sidebar.header('XRD anomolous scattering')
status_text = st.sidebar.empty()
#eta_input=st.sidebar.number_input('Electrode cover fraction (normalized)',0.0,1.0,0.3) #min, max, default
#t_input=st.sidebar.number_input('Thickness of layer (normalized)',0.0,10.0,0.2) #min, max, default

#st.sidebar.markdown('Beam')
#st.sidebar.markdown('h,k,l')
energy_kev = st.sidebar.number_input('Energy [keV]',0.0,1000000.0,17.0) #min, max, default

#st.sidebar.markdown('Geometry')
#sample_height = st.sidebar.number_input('Sample height [µm]',0.001,1000000.0,5000.0) #min, max, default
#sample_width = st.sidebar.number_input('Sample width [µm]',0.001,1000000.0,5000.0) #min, max, default
#sample_thickness = st.sidebar.number_input('Sample thickness [µm]',0.0,1000000.0,150.0) #min, max, default


import sys, os
import matplotlib.pyplot as plt
import numpy as np

uploaded_file = st.sidebar.file_uploader("Upload a .cif file")


known_crystal_structures = {
    'Diamond':'http://www.crystallography.net/cod/9012290.cif',
    'Aluminium':'http://www.crystallography.net/cod/9008460.cif',
    'BaTiO3':'http://www.crystallography.net/cod/1507756.cif',
    'KNbO3':'http://www.crystallography.net/cod/2310011.cif',
    'LiNbO3':'http://www.crystallography.net/cod/1541936.cif',
    'W':'http://www.crystallography.net/cod/9006486.cif',
    'Pb':'http://www.crystallography.net/cod/1531228.cif',
}
import os.path
import requests
keys = ['Upload']+list(known_crystal_structures.keys())
crystal = st.sidebar.selectbox('Or select a crystal structure',
        keys)

if not os.path.exists('cif_files'):
    os.makedirs('cif_files')
if uploaded_file is not None or crystal != 'Upload':
    if crystal == 'Upload':
        with open('cif_files/'+uploaded_file.name,"wb") as f:
            f.write(uploaded_file.getbuffer())
        cif_file = 'cif_files/'+uploaded_file.name #'KNO 2310011.cif'
        url = cif_file.split('/')[-1]
    else:
        url = known_crystal_structures[crystal]
        cif_file = 'cif_files/'+url.split('/')[-1]
        if not os.path.isfile(cif_file):
            r = requests.get(url)
            open(cif_file , 'wb').write(r.content)
    st.sidebar.write(f'Using {url}')



    import Dans_Diffraction
    from structure_factor import structure_factor
    lmbd = 4.135667696*10**-15 *299792458 / energy_kev /1000*10**6 #7.2932e-5 µm #7.2932e-11 m = 17 keV
    xtl = structure_factor.setup_dans_diffraction(cif_file, lmbd, ionize_fn = None)


    #xtl = Dans_Diffraction.Crystal(cif_file)
    #xtl.generate_lattice()
    #xtl.Scatter.setup_scatter(type='xray', energy_kev=energy_kev)


    material_str = ''.join(xtl.cif['_chemical_formula_sum'].replace('2','$_2$').replace('3','$_3$').split(' '))
    st.write(f'**{material_str}**  \n Density: {xtl.Properties.density():.3f} gm/cm'+r'$^3$')
    cell_par = str(xtl.Cell)
    cell_par = cell_par.replace('\n','').replace('\r','')
    cell_par = cell_par.split('A = ')
    #cell_par = [','.join(cps[0:3]), ','.join(cps[3:])]
    st.write(cell_par[0].replace('A','Å').replace('a','*a*').replace('b','*b*').replace('c','*c*'))
    st.write('𝛼 = '+cell_par[1].split('Volume')[0].replace('A','𝛼').replace('B','𝛽').replace('G','𝛾'))
    # pandas table
    xtl.Scatter._scattering_max_twotheta = 45
    reflections = xtl.Scatter.print_all_reflections().split('\n')
    #st.write(reflections[0])
    #st.write(reflections[2]+'  d_spacing [Å]')
    li = []
    for reflection in reflections[3:-2]:
        spli = reflection.split()
        #lambda = 2dsin(theta)

        d_spacing = lmbd*10**4/(2*np.sin(np.pi/360*float(spli[-2])))
        spl = reflection.split()
        li.append([' '.join(spl[:-2]), float(spl[-2]), float(spl[-1]),float(f'{d_spacing:13.3f}')])
        #st.write(reflection+f'  {d_spacing:13.3f}')


    import pandas as pd
    spl = reflections[2].split()
    df = pd.DataFrame(li, columns = [' '.join(spl[:-2]), '2𝜃', spl[-1]+' [%]', 'd-spacing [Å]'])

    feature_name = spl[-1]+' [%]'
    df[feature_name] = 100*df[feature_name] / df[feature_name].max()
    st.write(df.style.format(
                formatter={spl[-1]+' [%]': "{:.1f}",
                           '2𝜃': "{:.3f}",
                           'd spacing [Å]': "{:.3f}",
                          }
                ))


    # Q input
    hkl_str = st.sidebar.text_input('Q vector (h, k, l)', '('+li[0][0][1:-1]+')') #min, max, default
    hkl_spli = hkl_str.strip('()').split(',')
    h = int(hkl_spli[0])
    k = int(hkl_spli[1])
    l = int(hkl_spli[2])
    hkl = np.array([h,k,l])



    # ionization of ions:
    uvw,atom_types,label,occ,uiso,mxmymz = xtl.Structure.get()
    Q = xtl.Cell.calculateQ(hkl)[0]
    structure_factor.get_f0_WaasKirf('O',Q)
    ion_dropdown = []
    inonze_dict = {}
    for i,ion in enumerate(set(atom_types)):
        #st.write('ion')
        pos_keys = [ion+'2-',ion+'1-',ion+'5+',ion+'4+',ion+'3+',ion+'2+',ion+'1+',ion]
        keys = []
        for key in pos_keys:
            if key in structure_factor.f0_WaasKirf_dict.keys():
                keys.append(key)
        inonze_dict[ion] = keys[0]
        ion_dropdown.append( st.selectbox(ion+' ionization',keys, key = ion+'_'+str(i)) )


    for box in ion_dropdown:
        inonze_dict[ion] = box

    def ionize_fn(inp):
        if inp in inonze_dict.keys():
            return inonze_dict[inp]
        return inp

    ################# derivatives:
    a_scatter = xtl.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]
    atom_types_with_index = []
    for k, atom in enumerate(atom_types):
        if list(atom_types).count(atom) == 1:
            atom_types_with_index.append(atom)
        else:
            atom_types_with_index.append(atom+str(list(atom_types[0:k]).count(atom)))
    print_derivatives = st.checkbox('Print derivatives?')
    if print_derivatives:
        for k, _ in enumerate(a_scatter):
            absF = np.abs(np.sum(a_scatter))
            angF = np.angle(np.sum(a_scatter))
            absf_j = np.abs(a_scatter[k])
            angf_j = np.angle(a_scatter[k])
            Q = xtl.Cell.calculateQ(hkl)[0]
            dphidR_deg_per_AA = -absf_j/absF*np.cos(angf_j-angF)*Q*180/np.pi
            dabsF2dR_per_AA = -2*absf_j*absF*np.sin(-angf_j+angF)*Q/absF**2*100
            st.latex(r'\nabla_{R_{'+f'{atom_types_with_index[k]}'+'}} \phi_{F_{hkl}} = '+f'[{dphidR_deg_per_AA[1]*10**-2:.4f},{dphidR_deg_per_AA[1]*10**-2:.4f},{dphidR_deg_per_AA[2]*10**-2:.4f}]'+r'\ ^\circ\text{pm}^{-1}')
            st.latex(r'\nabla_{R_{'+f'{atom_types_with_index[k]}'+'}} |F_{hkl}|^2 = '+f'[{dabsF2dR_per_AA[1]*10**-2:.4f},{dabsF2dR_per_AA[1]*10**-2:.4f},{dabsF2dR_per_AA[2]*10**-2:.4f}]'+r'\ \%\text{pm}^{-1}')



    #################
    ################ move ions
    #################
    st.write('# Atomic posisions')
    move_options = ['-','Electric field','Manual' ]
    move_ions_selectbox = st.selectbox('Move ions?',move_options )

    if move_ions_selectbox == 'Electric field':
        default_weights_string = 'w'+str(atom_types_with_index[0])+'=1, '
        for atom in atom_types_with_index[1:]:
            default_weights_string += 'w'+str(atom)+'=0, '

        weights_string = st.text_input('', default_weights_string[:-2]) #min, max, default
        weights = np.zeros(len(atom_types_with_index))
        for i, s in enumerate(weights_string.split(',')):
            atom = s.split('=')[0].strip()
            index = atom_types_with_index.index(atom[1:])
            weights[index] = float( s.split('=')[1].strip())
        vareps_string = st.text_input('','𝛆ᵣ = [3600, 3600, 150]') #min, max, default
        vareps_temp = vareps_string.split('[')[1].split(']')[0].split(',')
        varpesilon= np.array([float(vareps_temp[0]),float(vareps_temp[1]),float(vareps_temp[2])])
        varpesilon_0 = 8.8541878128*10**-12
        Vuc = xtl.Cell.volume()*10**-30 #m^3
        charge_atoms = []
        for atom in atom_types:
            ionized = ionize_fn(atom).strip(atom).strip('+')
            if '-' in ionized:
                ionized = '-'+ionized[:-1]
            charge_atoms.append(float(ionized)*1.602176634*10**-19)
        charge_atoms = np.array(charge_atoms)

        Q = xtl.Cell.calculateQ(hkl)[0]
        Q_norm = 1/Q
        Q_norm[Q_norm==np.inf] = 0
        Q_norm = Q_norm/np.sqrt(np.sum(Q_norm**2))
        E_field_string = st.text_input('',f'E = [{Q_norm[0]:.3f},{Q_norm[1]:.3f},{Q_norm[2]:.3f}] kv/cm' ) #min, max, default
        E_field_temp = E_field_string.split('[')[1].split(']')[0].split(',')
        E_field = np.array([float(E_field_temp[0]),float(E_field_temp[1]),float(E_field_temp[2])])*1000/0.01 # V/m
        xtl_m = copy.deepcopy(xtl)
        for i, atom in enumerate(atom_types):
            if weights[i]>0:
                dr = weights[i]*(varpesilon-1)*varpesilon_0*Vuc/charge_atoms[i]*E_field*10**10
                uvw = np.array([xtl_m.Structure.u[i]*xtl.Cell.a+dr[0],
                                xtl_m.Structure.v[i]*xtl.Cell.b+dr[1],
                                xtl_m.Structure.w[i]*xtl.Cell.c+dr[2]])
                xtl_m.Structure.u[i] = uvw[0]/xtl.Cell.a
                xtl_m.Structure.v[i] = uvw[1]/xtl.Cell.b
                xtl_m.Structure.w[i] = uvw[2]/xtl.Cell.c
                if 0:
                    st.write(atom, charge_atoms[i])
                    st.write(dr)
                    absF = np.abs(np.sum(a_scatter))
                    angF = np.angle(np.sum(a_scatter))
                    absf_j = np.abs(a_scatter[i])
                    angf_j = np.angle(a_scatter[i])
                    dphidR_deg_per_AA = -absf_j/absF*np.cos(angf_j-angF)*Q*180/np.pi
                    dabsF2dR_per_AA = -2*absf_j*absF*np.sin(-angf_j+angF)*Q*100/absF/absF
                    st.write(f'dphidR = {np.sum(dphidR_deg_per_AA*dr):.2f}')
                    st.write(f'dabsF2dR = {np.sum(dabsF2dR_per_AA*dr):.2f}')

        a_scatter = xtl.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]
        a_scatter2 = xtl_m.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]

        st.write(f'Change in phase = {(np.angle(np.sum(a_scatter2))-np.angle(np.sum(a_scatter)))*180/np.pi:.4f} degrees')
        st.write(f'Change in |F|^2 = {(np.abs(np.sum(a_scatter2))**2-np.abs(np.sum(a_scatter))**2) /(np.abs(np.sum(a_scatter))**2)*100:.4f} %')

        if 0:
            Vuc =    3.9999*3.9999*4.0170*10**-30 #m
            qti = 4*1.602176634*10**-19 # C
            varpesilon = np.array([3600, 3600, 150])
            onekvpcm = 1*1000/0.01 # V/m
            E001 = np.array([0,0,1]) *onekvpcm
            E011 = np.array([0,1,1])/np.sqrt(2) *onekvpcm
            E111 = np.array([1,1,1])/np.sqrt(3) *onekvpcm
            delta_R_Ti_001 = (varpesilon-1)*varpesilon_0*Vuc/qti*E001*10**10
            delta_R_Ti_011 = (varpesilon-1)*varpesilon_0*Vuc/qti*E011*10**10
            delta_R_Ti_111 = (varpesilon-1)*varpesilon_0*Vuc/qti*E111*10**10

            st.write(f'delta_R_Ti_001 Å {delta_R_Ti_001}')
            st.write(f'delta_R_Ti_011 Å {delta_R_Ti_011}')
            st.write(f'delta_R_Ti_111 Å {delta_R_Ti_111}')
            dr = delta_R_Ti_111
        if 0:
            a_scatter = xtl.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]
            absF = np.abs(np.sum(a_scatter))
            angF = np.angle(np.sum(a_scatter))
            absf_j = np.abs(a_scatter[k])
            angf_j = np.angle(a_scatter[k])
            Q = xtl.Cell.calculateQ(hkl)[0]
            dphidR_deg_per_AA = -absf_j/absF*np.cos(angf_j-angF)*Q*180/np.pi
            dabsF2dR_per_AA = -2*absf_j*absF*np.sin(-angf_j+angF)*Q
            if 0:
                st.write(f'dphidR_deg_per_AA = {np.sum(np.sqrt(dphidR_deg_per_AA**2))}')
                st.write(f'{dphidR_deg_per_AA:.2f}')
                st.write(f'dabsF2dR_per_AA = {np.sum(np.sqrt(dabsF2dR_per_AA**2))}')
                st.write(f'{dabsF2dR_per_AA:.4f}')
            st.write(f'{np.sum(dphidR_deg_per_AA*dr):.4f} degrees')
            st.write(f'{np.sum(dabsF2dR_per_AA*dr)/(np.abs(np.sum(a_scatter))**2)*100:.4f} %')
            a_scatter = xtl.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]
            a_scatter2 = xtl_m.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]

            st.write(f'{(np.angle(np.sum(a_scatter2))-np.angle(np.sum(a_scatter)))*180/np.pi} degrees')
            st.write(f'{(np.abs(np.sum(a_scatter2))**2-np.abs(np.sum(a_scatter))**2) /(np.abs(np.sum(a_scatter))**2)*100} %')
            st.write(f'{np.abs(np.sum(a_scatter))**2/np.abs(np.sum(a_scatter))**2}')
            st.write(f'{np.abs(np.sum(a_scatter2))**2/np.abs(np.sum(a_scatter))**2}')

    elif move_ions_selectbox == 'Manual':
        xtl_m = copy.deepcopy(xtl)
        j = 0
        move_dropdowns = []
        atom_types = list(atom_types)
        atom_keys = ['-']+atom_types
        dr = np.array([0, 0, 0])

        while True:
            move_dropdowns.append( st.selectbox('Move ion [Å]',atom_keys, key = 'move_'+str(j)) )
            j+=1
            if move_dropdowns[-1] == '-':
                break
            else:
                k = atom_types.index(move_dropdowns[-1])
                #uvw = st.text_input('', 'u = '+str(xtl_m.Structure.u[k]*xtl.Cell.a)+', v = '+str(xtl_m.Structure.v[k]*xtl.Cell.b)+', w = '+str(xtl_m.Structure.w[k]*xtl.Cell.c),
                #    key = 'move_text_'+str(j)) #min, max, default

                uvw = st.text_input('', 'R = '+str(xtl_m.Structure.u[k]*xtl.Cell.a+dr[0])+', '+str(xtl_m.Structure.v[k]*xtl.Cell.b+dr[1])+', '+str(xtl_m.Structure.w[k]*xtl.Cell.c+dr[2]),
                    key = 'move_text_'+str(j)) #min, max, default
                xtl_m.Structure.u[k] = float(uvw.split('=')[1].split(',')[0])/xtl.Cell.a
                xtl_m.Structure.v[k] = float(uvw.split(',')[1])/xtl.Cell.b
                xtl_m.Structure.w[k] = float(uvw.split(',')[2])/xtl.Cell.c
                atom_keys.pop(k+1)
            a_scatter = xtl.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]
            a_scatter2 = xtl_m.Scatter.x_ray_anomalous(hkl, energy_kev = energy_kev, return_chi = False, per_atom = True, ionize_fn= ionize_fn)[0]

            st.write(f'{(np.angle(np.sum(a_scatter2))-np.angle(np.sum(a_scatter)))*180/np.pi} degrees')
            st.write(f'{(np.abs(np.sum(a_scatter2))**2-np.abs(np.sum(a_scatter))**2) /(np.abs(np.sum(a_scatter))**2)*100} %')

    if not 'xtl_m' in locals():
        fig, axes = plt.subplots(1,2,figsize=(8,4), dpi=100)
        axes = np.atleast_1d(axes)
        plot_scatter(xtl, axes[0], energy_kev= energy_kev)
        axes[0].set_title(f'Energy = {energy_kev:.2f} keV')
        axes[1].set_axis_off()
        st.write(fig)
        xtl_m = None
    else:
        fig, axes = plt.subplots(1,2,figsize=(8,4), dpi=100)
        axes = axes.ravel()
        plot_scatter(xtl, axes[0], energy_kev = float(energy_kev))
        axes[0].set_title(f'Energy = {energy_kev:.2f} keV')
        plot_scatter(xtl_m, axes[1], energy_kev= float(energy_kev))
        axes[1].set_title(f'Energy = {energy_kev:.2f} keV, moved ions')
        axes[1].set_ylabel('')
        st.write(fig)
    ################
    ################# select edge
    ################
    edges, energies = xtl.Properties.xray_edges()
    edge_keys = ['-']
    Q = xtl.Cell.calculateQ(hkl)[0]*10**4
    for edge, energy in zip(edges,energies):
        lmbd_e = 4.135667696*10**-15 *299792458 / energy /1000*10**6 #7.2932e-5 µm #7.2932e-11 m = 17 keV
        twotheta = 2*np.arcsin(np.sqrt(np.sum(Q**2))*lmbd_e/4/np.pi)*180/np.pi
        edge_keys.append(edge+f' {energy:.2f} keV, 2𝜃 = {twotheta:.2f}')
    edge_dropdown = st.selectbox('around resonant peak?',edge_keys)
    if not edge_dropdown == '-':
        i = edge_keys.index(edge_dropdown)-1
        E1 = energies[i]-0.01
        E2 = energies[i]+0.01
        fig, axes = plt.subplots(2,2,figsize=(8,8), dpi=100)
        axes = axes.ravel()
        axes = axes.ravel()
        plot_scatter(xtl, axes[0], energy_kev = float(E1))
        axes[0].set_title(f'Energy = {E1:.2f} keV')
        plot_scatter(xtl, axes[2], energy_kev= float(E2))
        axes[2].set_title(f'Energy = {E2:.2f} keV')
        axes[0].set_xlabel('')
        if not type(xtl_m) == type(None):
            plot_scatter(xtl_m, axes[1], energy_kev = float(E1))
            axes[1].set_title(f'Energy = {E1:.2f} keV, moved ions')
            plot_scatter(xtl_m, axes[3], energy_kev= float(E2))
            axes[3].set_title(f'Energy = {E2:.2f} keV, moved ions')
            axes[1].set_ylabel('')
            axes[1].set_xlabel('')
            axes[3].set_ylabel('')
        else:
            axes[1].set_axis_off()
            axes[3].set_axis_off()

        st.write(fig)
st.write("This toolbox is in beta testing, please independently verify all results, and let me know of any bugs :)")
st.write("tmara@dtu.dk")
