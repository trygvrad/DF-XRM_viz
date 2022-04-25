import streamlit as st
import numpy as np
import fpdf
import plotly
import sys, os
import os.path
import requests
import Dans_Diffraction
import pandas as pd

import plotly.io as pio
from PIL import Image

import helpers
import draw_crystal_structures
import three_d_draw
import three_d_draw_object_classes

import importlib

importlib.reload(sys.modules['three_d_draw_object_classes'])
importlib.reload(sys.modules['three_d_draw'])

# Yes. Change it to reportview-container - > appview-container
st.markdown(
        f"""
<style>
    .appview-container .main .block-container{{
        max-width: {1000}px;
        padding-top: {10}rem;
        padding-right: {2}rem;
        padding-left: {2}rem;
        padding-bottom: {10}rem;
    }}
    .appview-container .main {{
        color: {"#000000"};
        background-color: {"#EAEAEA"};
    }}
</style>
""",
        unsafe_allow_html=True,
    )

if 'key' not in st.session_state:
    st.session_state['key'] = 'value'

#######################################################################
################################# sidebar #############################
#######################################################################
st.header('DF-XRM visualization')
status_text = st.sidebar.empty()

def update_kev():
    st.session_state.energy_kev = 4.135667696*10**-15 *299792458 / st.session_state.lmbd_AA /1000*10**10

energy_kev = st.sidebar.number_input('Energy [keV]',0.0,1000000.0,17.0, key = 'energy_kev',format="%.4f") #min, max, default
lmbd_AA = st.sidebar.number_input('Wavelength [Å]',0.0,1000000.0,4.135667696*10**-15 *299792458 / energy_kev /1000*10**10, key = 'lmbd_AA', on_change = update_kev, format="%.4f")

st.sidebar.markdown('Geometry')
sample_thickness = st.sidebar.number_input('Sample thickness [µm]',0.0,1000000.0,150.0) #min, max, default

sample_height = st.sidebar.number_input('Sample height [µm]',0.001,1000000.0,5000.0) #min, max, default
sample_width = st.sidebar.number_input('Sample width [µm]',0.001,1000000.0,5000.0) #min, max, default


#######################################################################
############################## file input #############################
#######################################################################

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

keys = ['Upload']+list(known_crystal_structures.keys())
crystal = st.selectbox('Or select a crystal structure',
        keys)

# get cif file
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

    #######################################################################
    ######################### collect properties ##########################
    #######################################################################

    # load standard parameters
    params = {}

    params['wavelength_in_um'] = 4.135667696*10**-15 *299792458 / energy_kev /1000*10**6 #7.2932e-5 #7.2932e-11 m = 17 keV

    shape_µm = np.array([sample_height,sample_width,sample_thickness])
    params['shape'] = shape_µm
    params['mid'] = np.array([0.5, 0.5]) # axis 0, axis 1

    ''' get scatter function based on dan's diffraction
        Dan's diffraction reads cif files and is used to calculate correct atomic positions for scattering
        Dan's diffraction does not include a function for anomalous scattering, but this has been made based on xrddb and is assigned separately
        '''

    xtl = Dans_Diffraction.Crystal(cif_file)
    xtl.generate_lattice()
    wavelength_in_um = float(params['wavelength_in_um'])
    wavelength = wavelength_in_um*10**-6 # in meters
    xtl.Scatter.setup_scatter(type='xray', energy_kev=energy_kev)

    #######################################################################
    ###################### show material properties #######################
    #######################################################################

    material_str = ''.join(xtl.cif['_chemical_formula_sum'].replace('2','$_2$').replace('3','$_3$').split(' '))
    st.write(f'**{material_str}**  \n Density: {xtl.Properties.density():.3f} gm/cm'+r'$^3$')
    cell_par = str(xtl.Cell)
    cell_par = cell_par.replace('\n','').replace('\r','')
    cell_par = cell_par.split('A = ')
    #cell_par = [','.join(cps[0:3]), ','.join(cps[3:])]
    st.write(cell_par[0].replace('A','Å').replace('a','*a*').replace('b','*b*').replace('c','*c*'))
    st.write('𝛼 = '+cell_par[1].split('Volume')[0].replace('A','𝛼').replace('B','𝛽').replace('G','𝛾'))

    #######################################################################
    ########################## attenuation plot ###########################
    #######################################################################
    fig2d = helpers.get_att_plot(xtl, energy_kev, sample_thickness)
    st.write(fig2d)

    #######################################################################
    ####################### hkl - 2 theta table ###########################
    #######################################################################

    xtl.Scatter._scattering_max_twotheta = 45
    reflections = xtl.Scatter.print_all_reflections().split('\n')
    #st.write(reflections[0])
    #st.write(reflections[2]+'  d_spacing [Å]')
    li = []
    for reflection in reflections[3:-2]:
        spli = reflection.split()
        #lambda = 2dsin(theta)
        d_spacing = params['wavelength_in_um']*10**4/(2*np.sin(np.pi/360*float(spli[-2])))
        spl = reflection.split()
        li.append([' '.join(spl[:-2]), float(spl[-2]), float(spl[-1]),float(f'{d_spacing:13.3f}')])
        #st.write(reflection+f'  {d_spacing:13.3f}')


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

    #######################################################################
    ################################ get Q ################################
    #######################################################################

    st.sidebar.markdown('Scattering vector')
    hkl_str = st.sidebar.text_input('Q vector (h, k, l)', '('+li[0][0][1:-1]+')') #min, max, default
    hkl_spli = hkl_str.strip('()').split(',')
    h = int(hkl_spli[0])
    k = int(hkl_spli[1])
    l = int(hkl_spli[2])
    params['hkl'] = np.array([h,k,l]) #'3, 3, 3'
    Q = xtl.Cell.calculateQ(params['hkl'])[0]

    #######################################################################
    ####### the rotation of the crystal with respect to the sample ########
    #######################################################################
    if 0:
        '''
        Internally this code operates on the sample geometry
        the crystal lattice is rotated with respect to the sample
        and finally both the sample and crystal lattice are rotated according to the incomming beam
        '''

    st.sidebar.markdown('Crystal-sample geometry')
    up_hkl_str = st.sidebar.text_input("Sample 'up' (h,k,l) or [u,v,w,]", '') #min, max, default
    front_hkl_str = st.sidebar.text_input('Exit surface ~(h,k,l) or ~[u,v,w,]', '') #min, max, default

    if not up_hkl_str == '':
        s = up_hkl_str.strip('[()]').split(',')
        up_hkl = np.array([float(s[0]),float(s[1]),float(s[2])])
        if '(' in up_hkl_str:
            up_dir = xtl.Cell.calculateQ(up_hkl)[0]
        else:
            up_dir = xtl.Cell.calculateR(up_hkl)[0]
    else:
        up_dir = Q
        st.write(f'Sample up direction = Q')

    if not front_hkl_str == '':
        s = front_hkl_str.strip('[()]').split(',')
        front_hkl = np.array([float(s[0]),float(s[1]),float(s[2])])
        if '(' in front_hkl_str:
            front_dir = xtl.Cell.calculateQ(front_hkl)[0]
        else:
            front_dir = xtl.Cell.calculateR(front_hkl)[0]
    else:
        front_dir = np.cross(Q,np.array([0,0,1]))
        if np.sum(front_dir)<0.0001:
            front_dir = np.cross(Q,np.array([0,1,0]))
        front_index = xtl.Cell.indexQ(front_dir)[0]
        front_index/=np.sqrt(np.sum(front_index**2))
        if front_index[0]%1 == 0 and front_index[1]%1 == 0 and front_index[2]%1 == 0:
            st.write(f'Sample exit surface = ({front_index[0]:.0f},{front_index[1]:.0f},{front_index[2]:.0f})')
        else:
            st.write(f'Sample exit surface = ({front_index[0]:.2f},{front_index[1]:.2f},{front_index[2]:.2f})')


    if up_dir[0]>0:
        z_rot = -np.arctan2(up_dir[1],up_dir[0])*180/np.pi
    else:
        z_rot = -np.pi/2*np.sign(up_dir[1])*180/np.pi
    y_rot = -np.arctan2(up_dir[2],np.sqrt(up_dir[0]**2+up_dir[1]**2))*180/np.pi

    z_rot_crystal = z_rot
    y_rot_crystal = y_rot
    x_rot_crystal = 0

    def crystal_rotation_function(hkl):
        '''
        function that describes rotation in cartesian coordinates of crystal structure or reciprocal space
        performes an in-place rotation on a numpy array of length three
        '''
        helpers.rotate_z(hkl, z_rot_crystal*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        helpers.rotate_y(hkl, y_rot_crystal*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        helpers.rotate_x(hkl, x_rot_crystal*np.pi/180)# 45 degree rotation around the x-axis -> swap b and c
        return hkl

    front_dir_r = crystal_rotation_function(np.copy(front_dir))
    x_rot = np.arctan2(front_dir_r[1],front_dir_r[2])*180/np.pi
    x_rot_crystal = x_rot

    Q = crystal_rotation_function(Q)*10**4 # in inverse um
    params['Q_vector'] =  Q

    #######################################################################
    ######## Get the rotation of the sample wrt. the instrument ###########
    #######################################################################
    if 0:
        '''
        The user may here specify a rotation of the sample around the Q vector.
        This is implemented by specifying of k0 and kh in the sample coordinate system.
        (i.e. k0 and kh lie in the plane of Q and "Beam exit direction")

        in most cases the "Beam exit direction" will be the same as the "Sample exit surface"
        '''
    lmbd = float(params['wavelength_in_um'])
    k_abs = 2*np.pi/lmbd
    Q = np.linalg.norm(params['Q_vector'])
    st.sidebar.markdown('Beam-crystal geometry')
    beam_exit_hkl_str = st.sidebar.text_input('Beam exit direction ~(h,k,l) or ~[u,v,w,]', '') #min, max, default
    st.sidebar.write("For the 'exit surface' and 'sampe up' vectors, only the component orthogonal to Q is perserved")

    if not beam_exit_hkl_str == '':
        s = beam_exit_hkl_str.strip('[()]').split(',')
        beam_exit_hkl = np.array([float(s[0]),float(s[1]),float(s[2])])
        if '(' in beam_exit_hkl_str:
            beam_exit_dir = xtl.Cell.calculateQ(beam_exit_hkl)[0]
        else:
            beam_exit_dir = xtl.Cell.calculateR(beam_exit_hkl)[0]
    else:
        beam_exit_dir = front_dir
        exit_index = xtl.Cell.indexQ(beam_exit_dir)[0]
        exit_index/=np.sqrt(np.sum(exit_index**2))
        if exit_index[0]%1 == 0 and exit_index[1]%1 == 0 and exit_index[2]%1 == 0:
            st.write(f'Beam exit direction = ({exit_index[0]:.0f},{exit_index[1]:.0f},{exit_index[2]:.0f})')
        else:
            st.write(f'Beam exit direction = ({exit_index[0]:.2f},{exit_index[1]:.2f},{exit_index[2]:.2f})')
    beam_exit_dir = crystal_rotation_function(beam_exit_dir) # transformation for crystal coordinate system to sample coordinate system
    plane_normal = np.cross(params['Q_vector'],beam_exit_dir)
    k_dir_orthogonal_to_Q = np.cross(plane_normal,params['Q_vector'])
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
    params['k0_vector'] = -0.5*params['Q_vector'] + x*k_dir_orthogonal_to_Q
    params['kh_vector'] = 0.5*params['Q_vector'] + x*k_dir_orthogonal_to_Q
    beam_plane_normal = np.cross(params['k0_vector'],np.cross(params['k0_vector'],params['Q_vector']))
    twotheta = 2*np.arcsin(np.sqrt(np.sum(params['Q_vector']**2))*params['wavelength_in_um']/4/np.pi)*180/np.pi
    st.write(f'Instrument alignment on [{"".join(hkl_spli)}]:  \n 2𝜃 = {twotheta:.3f}')

    #######################################################################
    ###################### make 3d visualization ##########################
    #######################################################################
    scale = 1.0
    extend_beam = [25000.0,25000.0]
    extend_imaging_system = (5000.0,7500.0)
    legend_pos_shift = [2,0,-4]
    crystal_scale = 600.0
    crystal_axes_shift = [-5,-3,-5]
    crysta_structure_position = [-7000, -1000, 6000]
    cage_list = ['Ti', 'Nb']
    oxygen_cage_radius = 2.5
    make_bonds = ['C','C']
    max_bond_length = 2.5
    min_bond_length = 0
    show_text = True
    lens_scale = 1.0
    params['transverse_width'] = 100.0
    params['beam_thickness'] = 10.0

    vo = st.sidebar.checkbox("Visualization options")
    if vo:
        def parse_vec(description, default_string):
            s = st.sidebar.text_input(description, default_string)
            return np.array([float(v) for v in s.split(',')])
        scale = st.sidebar.number_input('Sample scale',0.00001,1000000.0,scale) #min, max, default
        extend_beam = parse_vec('Extend beam front/back µm', str(extend_beam).strip('[]'))*scale
        extend_imaging_system = np.array([5000,st.sidebar.number_input('sample-lens distance µm',0.00001,1000000.0,extend_imaging_system[1])])*scale
        lens_scale = st.sidebar.number_input('lens scale (scaled to sample)',0.00001,1000000.0,lens_scale) #min, max, default

        params['transverse_width'] = st.sidebar.number_input('beam width, µm',0.00001,1000000.0,params['transverse_width']) #min, max, default
        params['beam_thickness'] = st.sidebar.number_input('beam thickness, µm',0.00001,1000000.0,params['beam_thickness']) #min, max, default
        crystal_scale = st.sidebar.number_input('crystal structure scale',0.00001,1000000.0,crystal_scale)
        crysta_structure_position = parse_vec('crystal structure position', str(crysta_structure_position).strip('[]'))
        legend_pos_shift = parse_vec('crystal legend position', str(legend_pos_shift).strip('[]'))
        crystal_axes_shift = parse_vec('crystal axes position', str(crystal_axes_shift).strip('[]'))
        cage_list = st.sidebar.text_input('oxygen cages around', ', '.join(cage_list)).split(',')
        oxygen_cage_radius = st.sidebar.number_input('oxygen cage radius',0.00001,1000000.0,oxygen_cage_radius) #min, max, default
        make_bonds = st.sidebar.text_input('make bonds from/to', ', '.join(make_bonds)).split(',')
        make_bonds = [b.strip() for b in make_bonds]
        bond_length = parse_vec('bond length', str([min_bond_length, max_bond_length]).strip('[]'))
        min_bond_length = bond_length[0]
        max_bond_length = bond_length[1]
        show_text = st.sidebar.checkbox("show text", show_text)

    def make_crystal_structure():
            return draw_crystal_structures.add_crystal_structure( cif_file, scale = crystal_scale,
                        rotation_function = crystal_rotation_function,
                        legend_pos_shift = legend_pos_shift,
                        axes_shift = crystal_axes_shift,
                        max_bond_length = max_bond_length, min_bond_length = min_bond_length,
                        linewidth = 3,
                        bounding_box_facecolor = [0.5,0.6,0.7,0],
                        cage_line_color = [0.8,0.8,0.8,1],
                        linecolor = [0.5,0.5,0.5,1],
                        show_text = show_text,
                        cage_list = cage_list,
                        make_bonds = make_bonds
                        )

    drawing, fig, ax, outstr = three_d_draw.make_3d_perspective(params,
                export_blender_filepath = 'DF_XRM_vis_to_blender.pickled',
                crystal_rotation_function = crystal_rotation_function,
                scale = scale,
                tilt_to_vector = beam_plane_normal,
                atom_list = make_crystal_structure(),
                crysta_structure_position = crysta_structure_position,
                extend_beam = extend_beam,
                beam_opacity= 1.0,
                extend_imaging_system = extend_imaging_system,
                imaging_system_lw = 0,
                lens_scale = lens_scale*0.2,
                bounding_box_facecolor = [0,0,0,0.3],
                draw_stage = True,
                draw_axes = True,
                draw_scattering_arrows = True,
                show_text = show_text,
                )

    for i in range(len(outstr)):
        st.write(outstr[i].replace('phi','𝜑').replace('chi','𝜒').replace('omega','𝜔'))

    st.plotly_chart(fig)
    #######################################################################
    ########################### print to pdf ##############################
    #######################################################################

    if 1:
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
        txt.append(f'2theta = {twotheta:.3f}')

        if not up_hkl_str == '':
            txt.append(f"Sample 'up' hkl {up_hkl_str}")
        else:
            txt.append(f'Sample up direction = Q')
        if not front_hkl_str == '':
            txt.append(f'Exit surface {front_hkl_str}')
        else:
            txt.append(f'Sample exit surface = ({front_index[0]:.2f},{front_index[1]:.2f},{front_index[2]:.2f})')
        for i in range(len(outstr)):
            txt.append(outstr[i])
        if not beam_exit_hkl_str == '':
            txt.append(f'Beam exit direction = {beam_exit_hkl_str}')
        else:
            txt.append(f'Beam exit direction = ({exit_index[0]:.2f},{exit_index[1]:.2f},{exit_index[2]:.2f})')

        txt = '\n'.join(txt)
        txt = txt.encode('latin-1',  'ignore').decode('latin-1')
        pdf.set_xy(20, 20)
        pdf.multi_cell(0, 4, txt, 0, 1,'L', False)

        fig.write_image("3dfig.png")


        scope = pio.kaleido.scope
        scope._shutdown_kaleido()

        fig2d.savefig("2dfig.png")
        pdf.image('2dfig.png', w=200)

        test_image = "3dfig.png"
        original = Image.open(test_image)
        width, height = original.size   # Get dimensions
        left = width/6
        top = height/4
        right = 5 * width/6
        bottom = 3 * height/4
        cropped_example = original.crop((left, top, right, bottom))
        cropped_example.save('3dfig_crop.png')

        pdf.image('3dfig_crop.png', w=200)


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

    st.markdown(get_binary_file_downloader_html('DF_XRM_vis_to_blender.pickled', 'pickled 3d structure for blender'), unsafe_allow_html=True)
    st.markdown(get_binary_file_downloader_html('blender_import_script.py', 'blender import script'), unsafe_allow_html=True)

    #######################################################################
    ######################## export blender script ########################
    #######################################################################
    #if 1:


st.write("This is version 1.0.0. Souce code available at https://github.com/trygvrad/DF-XRM_viz under MIT lisence.")
st.write("please report bugs to tmara@dtu.dk, or contribute directly at github.")
