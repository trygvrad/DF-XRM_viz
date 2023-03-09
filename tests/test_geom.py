import numpy as np
import pandas as pd
import sys, os
import requests
import Dans_Diffraction
import plotly
import toml
import fpdf
import plotly.io as pio
import PIL.Image

from libs import absorption_plots
from libs import draw_crystal_structures
from libs import three_d_draw
from libs import optics_geometry
from libs import orientation_helpers

def get_cif_file(url):
    if not os.path.exists('cif_files'):
        os.makedirs('cif_files')
    cif_file = 'cif_files/'+url.split('/')[-1]
    if not os.path.isfile(cif_file):
        r = requests.get(url)
        open(cif_file , 'wb').write(r.content)
    return cif_file




def test_geometry():
    # input
    energy_kev = 17
    sample_thickness = 150
    params = {}
    params['wavelength_in_um'] = 4.135667696*10**-15 *299792458 / energy_kev /1000*10**6 #7.2932e-5 #7.2932e-11 m = 17 keV
    params['shape'] = np.array([5000,5000,sample_thickness])
    url = 'http://www.crystallography.net/cod/1533978.cif'
     
    crystal = 'MnO3Y'
    
    hkl_str = '(1,0,1)'
    up_hkl_str = '(0,0,1)'
    front_hkl_str = '(0,1,0)'
    
    is_beam_norm = False# "Minimize distance beam travels through sample" 
    beam_exit_hkl_str = '(0,1,0)' # only used if is_beam_norm == False
    
    lens_file_path = f'assets/lens_files/id06.lens'
    mainx = "d_tot' = 5000" # sample-detector distance, "d_tot or d_tot' [mm]"

    #############
    # answers for tests
    Q_answer = np.array([5490.7589723,  10219.38634624,  5900.16545805])
    two_theta_answer = 8.664149507727018

    z_rot_answer = 0
    y_rot_answer = -90
    x_rot_answer = 120
    ##############
    
    
    # setup 
    cif_file = get_cif_file(url)
    xtl = Dans_Diffraction.Crystal(cif_file)
    xtl.generate_lattice()
    xtl.Scatter.setup_scatter(scattering_type='xray', energy_kev=energy_kev)
    
    material_str = ''.join(xtl.cif['_chemical_formula_sum'].replace('2','$_2$').replace('3','$_3$').split(' '))
    
    cell_par = str(xtl.Cell)
    cell_par = cell_par.replace('\n','').replace('\r','')
    cell_par = cell_par.split('A = ')
    # make absorption fig    
    fig2d = absorption_plots.get_att_plot(xtl, energy_kev, sample_thickness)
    fig2d.savefig('2dfig.png')
    
    # scattering table 
    xtl.Scatter._scattering_max_twotheta = 45
    reflections = xtl.Scatter.print_all_reflections().split('\n')
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
    pd.set_option('display.max_rows', 300)
    tab = str(df)
    
    # orientation
    hkl_spli = hkl_str.strip('()').split(',')
    h = int(hkl_spli[0])
    k = int(hkl_spli[1])
    l = int(hkl_spli[2])
    params['hkl'] = np.array([h,k,l]) #'3, 3, 3'
    Q = xtl.Cell.calculateQ(params['hkl'])[0]
    
    up_dir = orientation_helpers.up_dir_from_string(xtl, up_hkl_str, Q)
    front_dir = orientation_helpers.front_dir_from_string(xtl, front_hkl_str, Q)
    
    z_rot = -np.arctan2(up_dir[1],up_dir[0])*180/np.pi
    y_rot = -np.arctan2(up_dir[2],np.sqrt(up_dir[0]**2+up_dir[1]**2))*180/np.pi

    
    x_rot = 0

    cr = orientation_helpers.crystal_rotation(z_rot, y_rot, x_rot)
    front_dir_r = cr.rotate(np.copy(front_dir))
    
    x_rot = np.arctan2(front_dir_r[1],front_dir_r[2])*180/np.pi
    
    cr = orientation_helpers.crystal_rotation(z_rot, y_rot, x_rot)
    crystal_rotation_function = cr.rotate

    
    print(f'z_rot, y_rot, x_rot = {z_rot}, {y_rot}, {x_rot}')
    assert np.allclose(np.array([z_rot, y_rot, x_rot]),\
           np.array([z_rot_answer, y_rot_answer, x_rot_answer])),  \
          f'np.array([z_rot, y_rot, x_rot]) vector not as expected '+\
          f'{np.array([z_rot, y_rot, x_rot])} expected'+\
          f' {np.array([z_rot_answer, y_rot_answer, x_rot_answer])}' 

    Q = crystal_rotation_function(Q)*10**4 # in inverse um
    params['Q_vector'] =  Q
    
    print(f'Q: {Q}')
    assert np.allclose(Q, Q_answer), f'Q vector not as expected {Q} expected {Q_answer}' 

    
    k_abs = 2*np.pi/params['wavelength_in_um']
    Q = np.linalg.norm(params['Q_vector'])
    
    '''is_beam_norm = st.sidebar.checkbox("Minimize distance beam travels through sample", True)
    beam_exit_hkl_str = ''
    if not is_beam_norm:
        beam_exit_hkl_str = st.sidebar.text_input('Beam exit direction ~(h,k,l) or ~[u,v,w,]', '') #min, max, default
    '''
    if not is_beam_norm and not beam_exit_hkl_str == '':
        beam_exit_dir = orientation_helpers.direction_from_string(xtl, beam_exit_hkl_str)
    else:
        beam_exit_dir = front_dir
        exit_index = xtl.Cell.indexQ(beam_exit_dir)[0]
        exit_index/=np.sqrt(np.sum(exit_index**2))
    beam_exit_dir = crystal_rotation_function(beam_exit_dir) # transformation for crystal coordinate system to sample coordinate system
    plane_normal = np.cross(params['Q_vector'],beam_exit_dir)
    k_dir_orthogonal_to_Q = np.cross(plane_normal,params['Q_vector'])
    k_ort_l = np.linalg.norm(k_dir_orthogonal_to_Q)

    
    x = np.sqrt(  ( k_abs**2 - 0.25*Q**2 ) / ( k_ort_l**2) )
    params['k0_vector'] = -0.5*params['Q_vector'] + x*k_dir_orthogonal_to_Q
    params['kh_vector'] = 0.5*params['Q_vector'] + x*k_dir_orthogonal_to_Q
    two_theta = 2*np.arcsin(np.sqrt(np.sum(params['Q_vector']**2))*params['wavelength_in_um']/4/np.pi)*180/np.pi
    print(f'two_theta: {two_theta}')
    assert two_theta == two_theta_answer, f'two_theta not as expected {two_theta} expected {two_theta_answers}' 

    # make 3d fig  -  default options
    scale = 1.0
    extend_beam = [25000.0,25000.0]
    extend_scattered_beam = 7500.0
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
    magnification = 3.0

    extend_scattered_beam *= 3 # when showing lens


    crystal_structure = draw_crystal_structures.add_crystal_structure( cif_file, scale = crystal_scale,
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
                make_bonds = make_bonds,
                #atom_legend_step = Q/np.sqrt(np.sum(Q**2)),
                )

    drawing, fig, ax, outstr = three_d_draw.make_3d_perspective(params,
                export_blender_filepath = 'DF_XRM_vis_to_blender.pickled',
                crystal_rotation_function = crystal_rotation_function,
                scale = scale,
                atom_list = crystal_structure,
                crysta_structure_position = crysta_structure_position,
                extend_beam = extend_beam,
                extend_scattered_beam = extend_scattered_beam,
                lens_scale = lens_scale*0.2,
                bounding_box_facecolor = [0,0,0,0.3],
                draw_stage = True,
                draw_axes = True,
                draw_scattering_arrows = True,
                show_text = show_text,
                magnification = magnification,
                )
    
    # make lens fig
    if "d_tot'" in mainx:
        d_tot = float(mainx.split("=")[1])/np.cos(two_theta*np.pi/180)
    else:
        d_tot = float(mainx.split("=")[1])
    with open(lens_file_path) as f:
        lens_file_contents = f.read()
    parsed_lens_file = toml.loads(lens_file_contents)
    fig_optics, ax = optics_geometry.make_optics_geometry_plot(\
                    energy_kev, two_theta, d_tot, parsed_lens_file)
    fig_optics.savefig('fig_optics.png')
    
    
    # make pdf 
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
    txt.append('A = ' +cell_par[1].replace('Volume',' Volume'))
    txt.append(f'Q vector hkl {hkl_str}')
    txt.append(f'2theta = {two_theta:.3f}')

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
    original = PIL.Image.open(test_image)
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
    if True: #has_lens_fig:
        pdf.add_page()
        pdf.set_xy(20, 20)
        pdf.multi_cell(0, 4, txt, 0, 1,'L', False)
        fig_optics.savefig("fig_optics.png")
        pdf.image('fig_optics.png', w=200)

    pdf.output('test.pdf')

