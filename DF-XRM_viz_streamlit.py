import streamlit as st
import time
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import plotly

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

st.header('3D-XRM visualization')
status_text = st.sidebar.empty()
#eta_input=st.sidebar.number_input('Electrode cover fraction (normalized)',0.0,1.0,0.3) #min, max, default
#t_input=st.sidebar.number_input('Thickness of layer (normalized)',0.0,10.0,0.2) #min, max, default

#st.sidebar.markdown('Beam')
#st.sidebar.markdown('h,k,l')
energy_kev = st.sidebar.number_input('Energy [keV]',0.0,1000000.0,17.0) #min, max, default

st.sidebar.markdown('Geometry')
sample_height = st.sidebar.number_input('Sample height [µm]',0.001,1000000.0,5000.0) #min, max, default
sample_width = st.sidebar.number_input('Sample width [µm]',0.001,1000000.0,5000.0) #min, max, default
sample_thickness = st.sidebar.number_input('Sample thickness [µm]',0.0,1000000.0,150.0) #min, max, default


if 0:
    from PIL import Image
    image = Image.open('streamlit IDE.png')

    st.image(image, caption='Notation used to describe the geometry of interdigitated electrodes on a thin film',
          use_column_width=True)
    st.subheader('Capacitance of the structure described by the input parameters:')

    st.write( C_tot, 'F')

import sys, os
import matplotlib.pyplot as plt
import numpy as np

uploaded_file = st.file_uploader("Upload a .cif file")


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
crystal = st.selectbox('Or select a crystal structure',
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
    st.write(f'Using {url}')

    # load standard parameters
    params = {}
    params['BEAM'] = {}
    params['GEOMETRY'] ={}
    params['BEAM']['transverse_width'] = 100

    # energy_kev = 4.135667696*10**-15 *299792458 / wavelength /1000 # h[eV⋅s]*c[m⋅s]/wavelength[m]
    #wavelength =
    params['BEAM']['lmbd'] = 4.135667696*10**-15 *299792458 / energy_kev /1000*10**6 #7.2932e-5 #7.2932e-11 m = 17 keV
    #print("params['BEAM']['lmbd']",params['BEAM']['lmbd'])
    params['GEOMETRY']['shape'] = np.array([300, 300, 10]) #[16000, 1600, 2000]
    shape_µm = np.array([sample_height,sample_width,sample_thickness])
    params['GEOMETRY']['del_x'] = shape_µm/params['GEOMETRY']['shape']#'0.05, 0.05, 0.05'
    params['BEAM']['mid'] = np.array([0.5, 0.5]) # axis 0, axis 1



    if 0:
        ''' get scatter function based on dan's diffraction
            Dan's diffraction reads cif files and is used to calculate correct atomic positions for scattering
            Dan's diffraction does not include a function for anomalous scattering, but this has been made based on xrddb and is assigned separately
            '''

    import Dans_Diffraction


    xtl = Dans_Diffraction.Crystal(cif_file)
    xtl.generate_lattice()
    wavelength_in_um = float(params['BEAM']['lmbd'])
    wavelength = wavelength_in_um*10**-6 # in meters
    xtl.Scatter.setup_scatter(type='xray', energy_kev=energy_kev)


    ########################## crystal rotation functions #####################
    def rotate(x,z,phi):# 2D rotation function
        '''
        generic rotation funciton for generating the crystal rotation function
        '''
        return x*np.cos(phi)-z*np.sin(phi),z*np.cos(phi)+x*np.sin(phi)
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
    #f1f2_BrennanCowan = pkg_resources.read_text(databases, 'f1f2_BrennanCowan.dat')


    material_str = ''.join(xtl.cif['_chemical_formula_sum'].replace('2','$_2$').replace('3','$_3$').split(' '))
    st.write(f'**{material_str}**  \n Density: {xtl.Properties.density():.3f} gm/cm'+r'$^3$')
    cell_par = str(xtl.Cell)
    cell_par = cell_par.replace('\n','').replace('\r','')
    cell_par = cell_par.split('A = ')
    #cell_par = [','.join(cps[0:3]), ','.join(cps[3:])]
    st.write(cell_par[0].replace('A','Å').replace('a','*a*').replace('b','*b*').replace('c','*c*'))
    st.write('𝛼 = '+cell_par[1].split('Volume')[0].replace('A','𝛼').replace('B','𝛽').replace('G','𝛾'))


    fig, axes = plt.subplots(1,2,figsize=(12,4), dpi=100)
    fig2d = fig
    import xraydb
    def att_length(xtl,energy_kev):
        atom_type = xtl.Structure.type
        occ = xtl.Structure.occupancy
        natoms = np.sum(occ)
        vol = xtl.Properties.volume()*10**6/10**30 #Å^3
        weight = Dans_Diffraction.fc.atom_properties(atom_type, 'Weight')/(6.02214076*10**23)
        #st.write(weight)
        mu =  occ[0]*weight[0]*xraydb.mu_elam(atom_type[0], energy_kev*1000)
        #mu +=  occ[0]*weight[0]*xraydb.incoherent_cross_section_elam(atom_type[0], energy_kev*1000)
        # Values returned are in units of cm^2/gr # here mutiplied by weight
        for i in range(1,len(occ)):
            mu +=  occ[i]*weight[i]*xraydb.mu_elam(atom_type[i], energy_kev*1000)#photoabsorption_crosssection(elements, energy_kev)
            #mu +=  occ[0]*weight[0]*xraydb.incoherent_cross_section_elam(atom_type[0], energy_kev*1000)
        #st.write(vol)
        absorption = vol/mu
        return absorption*10**4 # cm to µm

    ax = axes[0]
    energies = np.arange(1,30,0.1)
    #absorption_invum = xtl.Properties.absorption(energy_kev=energies) #Returns the sample absorption coefficient in um^-1 at the requested energy in keV
    #attenuation_um = xtl.Properties.xray_attenuation_length(energy_kev=energies) # the *10 == bug in dan's diffraction?
    attenuation_um2 = att_length(xtl, energies) # the *10 == bug in dan's diffraction?

    #ax.loglog(energies, attenuation_um, color = [0.8,0.2,1])
    ax.loglog(energies, attenuation_um2, color = [0.8,0.2,0])
    ax.set_xlabel('energy [keV]')
    ax.set_ylabel(r'Attenuation length [µm]')
    #ax.set_xticks([1,2,3,4,5,6,7,8,9,10,20,30], minor = True)
    #ax.set_xticklabels([1,2,3,4,'',6,'',8,'',10,20,30], minor = True)
    ax.set_xticks([1,10])
    ax.set_xticklabels([1,10])
    y_lim = ax.get_ylim()
    ax.plot([energy_kev, energy_kev], y_lim, color = [0,0.3,1,0.5])
    ax.set_ylim(y_lim)
    ax.text(energy_kev,y_lim[0]*2, f'E = {energy_kev:.2f} keV', rotation = 'vertical', ha = 'right')

    att = att_length(xtl, xtl.Scatter._energy_kev)
    #att = xtl.Properties.xray_attenuation_length(energy_kev = xtl.Scatter._energy_kev)[0]
    ax.text(energy_kev,att, f'{att:.0f} µm', ha = 'right')



    ax = axes[1]
    thickness = np.arange(50,1000,1)
    #absorption_invum = xtl.Properties.absorption(energy_kev=energies) #Returns the sample absorption coefficient in um^-1 at the requested energy in keV
    #att = xtl.Properties.xray_attenuation_length(energy_kev= xtl.Scatter._energy_kev)[0]*10
    xray_transmission = np.exp(-thickness/att)
    ax.loglog(thickness, xray_transmission, color = [0.8,0.2,0])
    ax.set_xlabel('Thickness [µm]')
    ax.set_ylabel(r'Transmission [-]')
    ax.set_title(f'Transmission @ {np.round(xtl.Scatter._energy_kev):.0f} keV')
    #ax.set_xticks([1,2,3,4,5,6,7,8,9,10,20,30], minor = True)
    #ax.set_xticklabels([1,2,3,4,'',6,'',8,'',10,20,30], minor = True)
    #ax.set_xticks([1,10])
    #ax.set_xticklabels([1,10])
    y_lim = ax.get_ylim()
    if y_lim[1]<10*y_lim[0]:
        y_lim=[0.1,1.1]
    ax.plot([sample_thickness, sample_thickness], y_lim, color = [0,0.3,1,0.5])
    ax.set_ylim(y_lim)
    ax.text(sample_thickness,y_lim[0]*2, f'Thickness = {sample_thickness:.2f} µm', rotation = 'vertical', ha = 'right')
    xray_transmission_0 = np.exp(-sample_thickness/att)
    ax.text(sample_thickness,xray_transmission_0, f'{xray_transmission_0*100:.2f}%')




    st.write(fig)



    xtl.Scatter._scattering_max_twotheta = 45
    reflections = xtl.Scatter.print_all_reflections().split('\n')
    #st.write(reflections[0])
    #st.write(reflections[2]+'  d_spacing [Å]')
    li = []
    for reflection in reflections[3:-2]:
        spli = reflection.split()
        #lambda = 2dsin(theta)
        d_spacing = params['BEAM']['lmbd']*10**4/(2*np.sin(np.pi/360*float(spli[-2])))
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
    #a_series = pd.Series(li, index = df.columns)
    #df = df.append(a_series, ignore_index=True)



    hkl_str = st.text_input('Q vector h, k, l', li[0][0][1:-1]) #min, max, default
    hkl_spli = hkl_str.split(',')
    h = int(hkl_spli[0])
    k = int(hkl_spli[1])
    l = int(hkl_spli[2])
    params['GEOMETRY']['hkl'] = np.array([h,k,l]) #'3, 3, 3'


    hkl = params['GEOMETRY']['hkl']
    Q = xtl.Cell.calculateQ(hkl)[0]




    def crystal_rotation_function(hkl):
        '''
        function that describes rotation in cartesian coordinates of crystal structure or reciprocal space
        performes an in-place rotation on a numpy array of length three
        '''
        rotate_z(hkl, z_rot_crystal*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        rotate_y(hkl, y_rot_crystal*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        rotate_x(hkl, x_rot_crystal*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        return hkl

    if True:#st.checkbox('Custom unit cell rotation'):
        up_hkl_str = st.text_input("Sample 'up' h, k, l", hkl_str) #min, max, default
        s = up_hkl_str.split(',')
        up_hkl = np.array([float(s[0]),float(s[1]),float(s[2])])
        if np.sum(np.abs(np.cross(up_hkl, np.array([0,0,1])))) >0:
            front_hkl_str = st.text_input('Exit surface h, k, l', '0, 0 ,1') #min, max, default
        else:
            front_hkl_str = st.text_input('Exit surface h, k, l', '0, 1 ,0') #min, max, default
        s = front_hkl_str.split(',')
        front_hkl = np.array([float(s[0]),float(s[1]),float(s[2])])

        up_dir = xtl.Cell.calculateR(up_hkl)[0]
        front_dir = xtl.Cell.calculateR(front_hkl)[0]

        if up_dir[0]>0:
            z_rot = -np.arctan2(up_dir[1],up_dir[0])*180/np.pi
        else:
            z_rot = -np.pi/2*np.sign(up_dir[1])*180/np.pi
        y_rot = -np.arctan2(up_dir[2],np.sqrt(up_dir[0]**2+up_dir[1]**2))*180/np.pi

        z_rot_crystal = z_rot
        y_rot_crystal = y_rot
        x_rot_crystal = 0
        front_dir = crystal_rotation_function(front_dir)
        x_rot = np.arctan2(front_dir[1],front_dir[2])*180/np.pi
        x_rot_crystal = x_rot


        #z_rot_crystal = st.number_input('Rotate unit cell around Z',-360.0,360.0,z_rot) #min, max, default
        #y_rot_crystal = st.number_input('Then rotate around Y',-360.0,360.0,y_rot) #min, max, default
        #x_rot_crystal = st.number_input('Then rotate around X',-360.0,360.0,x_rot) #min, max, default



    Q = crystal_rotation_function(Q)*10**4 # in inverse um
    params['BEAM']['Q_vector'] =  Q

    import importlib
    import draw_crystal_structures
    importlib.reload(sys.modules['draw_crystal_structures'])
    def make_crystal_structure():
            return draw_crystal_structures.add_crystal_structure( cif_file, scale = 6.0,
                                            rotation_function = crystal_rotation_function,
                                            legend_pos_shift = [0,-50,-50],
                                            axes_shift = [-5,-3,-5],
                                            #make_bonds=['O','O '],
                                            max_bond_length = 1.5, min_bond_length = 0,
                                            linewidth = 3,
                                            bounding_box_facecolor = [0.5,0.6,0.7,0],
                                            cage_line_color = [0.8,0.8,0.8,1],
                                            linecolor = [0.5,0.5,0.5,1])
    # import the renderer
    import three_d_draw
    import three_d_draw_object_classes
    import importlib
    importlib.reload(sys.modules['three_d_draw_object_classes'])
    importlib.reload(sys.modules['three_d_draw'])
    # what we want is a blender backend. Unfortunately blender does not install easily as a module (without gui).
    # Therefore use mayavi for 3d rendering (slow) or matplotlib for pseudo-3d rendering(fast)
    backend = 'plotly' # default is 'matplotlib'
    #backend = 'matplotlib'
    three_d_draw.set_backend(backend)
    #dfxrm.three_d_draw.timing=True # <- enable this if you want to time the rendering


    def final_rotation_function(drawing):
        #drawing.rot_y(-10*np.pi/180) # z in coordinate system of drawing is out of plane
        #drawing.rot_x(10*np.pi/180) # z in coordinate system of drawing is out of plane
        return

    import types
    ex = types.SimpleNamespace()
    ex.params = params
    ex.Q_vector = params['BEAM']['Q_vector'] # in inverse µm
    ex.shape = params['GEOMETRY']['shape']
    ex.del_x = params['GEOMETRY']['del_x']
    if 1:
        if 0:
                '''
                Calculate self.k0_vector and self.kh_vector so that they lie in the plane defined by self.Q_vector and the [001] axis
                input: None
                Return: None
                '''
        lmbd = float(params['BEAM']['lmbd'])
        k_abs = 2*np.pi/lmbd
        Q = np.linalg.norm(ex.Q_vector)
        plane_normal = np.cross(ex.Q_vector,np.array([0,0,1]))
        k_dir_orthogonal_to_Q = np.cross(plane_normal,ex.Q_vector)
        k_ort_l = np.linalg.norm(k_dir_orthogonal_to_Q)
        if 0:
            '''
            We have:
                k0 = -0.5*Q + x*k_ort
                so that |k0| = k
                where Q and k_ort are orthogonal, so that
                0.25|Q|**2 + x**2*|k_ort|**2 = |k0|**2
                i.e:
                x = np.sqrt(  ( |k0|**2 - 0.25|Q|**2 ) / ( |k_ort|**2) )
            '''
        x = np.sqrt(  ( k_abs**2 - 0.25*Q**2 ) / ( k_ort_l**2) )
        ex.k0_vector = -0.5*ex.Q_vector + x*k_dir_orthogonal_to_Q
        ex.kh_vector = 0.5*ex.Q_vector + x*k_dir_orthogonal_to_Q

    beam_plan_normal = np.cross(ex.k0_vector,np.cross(ex.k0_vector,ex.Q_vector))

    twotheta = 2*np.arcsin(np.sqrt(np.sum(ex.Q_vector**2))*params['BEAM']['lmbd']/4/np.pi)*180/np.pi
    st.write(f'Instrument alignment on [{"".join(hkl_spli)}]:  \n 2𝜃 = {twotheta:.3f}')

    if 1:
        drawing, fig, ax, = three_d_draw.make_3d_perspective(ex, #volume_data = make_internals(step_size = 1),
                                                                scale = 0.1/8,
                                                                beam_step_size = 1, show_beam_at_end = False,
                                                                view_angle=20, fig_size = (12,6), factor = 350,
                                                                tilt_to_vector = beam_plan_normal,
                                                                extend_imaging_system = (10*150/sample_thickness,50*150/sample_thickness),
                                                                imaging_system_lw = 0,
                                                                extend_beam = [30*150/sample_thickness,30*150/sample_thickness],
                                                                draw_axes = True,
                                                                draw_curved_arrows = True,
                                                                beam_draw_mesh = False,
                                                                #draw_face_axes = True,
                                                                crystal_rotation_function = crystal_rotation_function,
                                                                beam_opacity= 1.0,
                                                                #export_blender_filepath = 'pickle_to_blender.pickled',
                                                                atom_list = make_crystal_structure(),
                                                                #lens_scale = 0.5,
                                                                final_rotation_function= final_rotation_function,
                                                                bounding_box_facecolor = [0,0,0,0.3],
                                                                   )#(600*4, 300*4))
        st.plotly_chart(fig)

        ############### print to pdf #################
        import fpdf
        pdf = fpdf.FPDF(orientation = 'P', unit = 'mm', format = 'A4')
        pdf.add_page()
        pdf.set_font('helvetica', '', 10)
        #pdf.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf', uni=True)
        #pdf.set_font('DejaVu', '', 14)
        pdf.set_text_color(0, 0, 0)
        txt =[]
        txt.append(crystal)
        txt.append(f'Using {url}')
        txt.append(f'{material_str}\nDensity: {xtl.Properties.density():.3f} gm/cm3')
        txt.append(cell_par[0])
        txt.append('A = ' +cell_par[1])
        txt.append(f'Q vector hkl {hkl_str}')
        txt.append(f"Sample 'up' hkl {up_hkl_str}")
        txt.append(f'Exit surface {front_hkl_str}')
        txt = '\n'.join(txt)
        txt = txt.encode('latin-1', 'ignore').decode('latin-1')
        pdf.set_xy(20, 20)
        pdf.multi_cell(0, 4, txt, 0, 1,'L', False)

        fig.write_image("3dfig.png")
        fig2d.savefig("2dfig.png")
        pdf.image('2dfig.png', w=200)
        pdf.image('3dfig.png', w=200)


        pdf.add_page()
        pd.set_option('display.max_rows', 300)
        tab = str(df)
        tab = tab.replace('𝜃','theta').encode('latin-1', 'ignore').decode('latin-1')
        pdf.multi_cell(0, 4, tab, 0, 1,'L', False)
        pdf.output('DF-XRM_vis.pdf')

        import base64
        def get_binary_file_downloader_html(bin_file, file_label='File'):
            with open(bin_file, 'rb') as f:
                data = f.read()
            bin_str = base64.b64encode(data).decode()
            href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">Download {file_label}</a>'
            return href
        st.markdown(get_binary_file_downloader_html('DF-XRM_vis.pdf', 'pdf'), unsafe_allow_html=True)
st.write("This toolbox is in beta testing, please independently verify all results, and let me know of any bugs :)")
st.write("tmara@dtu.dk")
