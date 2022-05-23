---

author:
- 'Trygve Magnus Ræder, tmara\@dtu.dk'
date: April 2022
title: 'DF-XRM\_viz'
---



# Summary

`DF-XRM_vis` is a visualization tool to aid in planning, execution, and
data analysis of dark-field X-ray microscopy experiments. The toolkit is built in python and has a streamlit interface, as well as a jupyterlab example. The highlight of the toolkit is generating three dimensional models of experimental geometry, rendered together with the unit cell of the material. By hosting the streamlit application on a publically avaiable adress, the toolkit becomes available also to researchers without programming experience.

# Statement of need

Dark-field X-ray microscopy (DF-XRM) is a recent development in the field of X-ray science\[ref\]. DF-XRM relies on Bragg diffraction, and
uses a lens to produce a real-space image on a detector. The use of a
focusing lens differentiates DF-XRM from other X-ray imaging techniques
such using Bragg diffraction (Bragg scanning probe X-ray microscopy, Bragg coherent diffraction imaging and Bragg X-ray ptychography). All the aforementioned
techniques fundamentally capture the same crystallographic information,
and have different advantages and disadvantages. DF-XRM has the distinct
advantage that the results are natively images with real space axes (as
opposed to images of reciprocal space). This advantage makes DF-XRM more
approachable of researchers that may have limited
familiarity with reciprocal space. However, as a technique based on
Bragg diffraction, there are a number of geometrical and
crystallographic constraints that must be fulfilled in order to capture
a DF-XRM image. Failure to apprehend or communicate the constraints has
been a liability in the practical application of DF-XRM.
`DF-XRM_vis` is designed to be an aid in assessing the geometrical
constraints in X-ray microscopy. `DF-XRM_vis` may also aid in data analysis, serve as a
teaching tool, or an aid when discussing DF-XRM more broadly.

A successful DF-XRM experiment requires that the following pitfalls be
avoided:

1. The sample absorption at the desired X-ray energy is too high and/or
   the sample is too thick.

2. The desired reflection and X-ray energy does not produce a
   diffraction angle suitable for the instrument.

3. The sample mount assumes a different crystallographic orientation
   than what the crystal actually has.

Avoiding the above pitfalls requires good communication of 3d geometry
between the researchers involved in the experiment. Traditionally, this
has been an (often challenging) exercise involving the use of a
whiteboard, but as more communication has shifted online, collaboration
on 3d geometry has become more difficult. I have personally observed
experiments loose days of beamtime to the issues listed above, as a result
of poor online communication in the lead-up to experiments.


While `DF-XRM_vis` is designed to address all the above points, we note that the novelty of `DF-XRM_vis` lies with adressing the 3rd point, which has not previously been adressed in any comparable tool. We will here describe the capabilities of the streamlit interface to `DF-XRM_vis` in order of increasing user input, which is the same as the order as the points listed above.

Sample absorption
---------------------------

Uploading a `.cif` file to `DF-XRM_vis` (Fig. [1](#fig:1){reference-type="ref" reference="fig:1"}a) `DF-XRM_vis` will plot the attenuation length with respect to the energy (Fig. [1](#fig:1){reference-type="ref" reference="fig:1"}d), and plot the transmission with respect to the sample thickness (Fig. [1](#fig:1){reference-type="ref" reference="fig:1"}e).

A handful of crystal-structures are available from the dropdown menu for easy access (as shown in the example).

![Attenuation and transmission in `DF-XRM_vis`. **a** Upload or selection of `.cif` file. **b** selection of X-ray energy. **c** selection of sample thickness. **d** plot of energy versus attenuation. **e** plot of thickness versus transmission.\label{fig:1}](df-xrm_1.png){ width=800 }



The selected X-ray energy can then be changed (Fig. [1](#fig:1){reference-type="ref" reference="fig:1"}b, in keV or Å), and
the thickness of the sample can be changed (Fig. [1](#fig:1){reference-type="ref" reference="fig:1"}c).

This gives an immediate assessment of the flux exiting the sample, thus determins the viability of utilizing the sample for DF-XRM.

Considering Q-vectors
---------------------

Based on the provided `.cif` file and chosen wavelength, `DF-XRM_vis` will provide a table of available reflections, their scattering angle and *d*-spacing.

![Example the table produced after a `.cif` file is chosen. The example describes KNbO$_3$ at 17 keV.\label{fig:2}](df-xrm_2.png){width=800}

A suitable reflection, as allowed by the instrument, can then be chosen. `DF-XRM_vis` will automatically select the brightest reflection as an initial guess.

Crystal-sample-instrument orientation
-------------------------------------

A 3d rendering of the experimental conditions will then be generated based on an educated guess of sample facets. The user may then change the sample facets using real or reciprocal lattice vectors.

The rendering shows the required alignment of the crystal at to fulfill the selected reflection. Users may freely rotate the rendering in three dimensional space. The rendering consists of 4 parts:

1. The sample, annotated with the dimensions along the differnt axes.

2. The crystal structure, shown in the orientation determined by the scattering condition and and in complience with the crystallographic axes of the sample.

3. The beam, lens, and scattered beam.

4. The goiniometer stage, annotated with the angles of rotation.


![Example visualization showing an x-cut LiNbO$_3$ wafer, aligned on the (1,-1,2) (brightest) reflection.\label{fig:3}](df-xrm_3.png){width=800}

Rendering the sample with realistic dimensions makes it easier to ensure that samples are positioned with the intended orientation during the experiment, and helps visualize past experiments. 

Showing the crystal structure together with the sample geometry ensures that the relationship between the crystal lattice and sample geometry is correct. The rendering intentionally uses a parallel projection. 
This makes it easier to compare the crystal axes to the axes of the sample. 

By showing the beam in the rendering, the intersection between the beam and sample can be understood. The detector image in DF-XRM is a projection of the intersection between the beam and sample, so that a rendering of this intersection is useful when inspecting DF-XRM images. 

The goiniometer position is also shown. For large angles of $\phi$, such as shown here, it is advised to build a custom sample holder to hold the sample prior to the experiment. For smaller angles, the values presented in `DF-XRM_vis` may be directly applied to the instrument as part of the alignment procedure.

Other features
-------------------------------------

A link is provided to download the current view as a pdf. This functionality is intended for inclusion in the logbook of an ongoing experiment.

Additionally, the three dimensional rendering may be downloaded as a pickled library and imported to blender by using a script provided.

Comparison to other tools
-------------------------

other tools for calculating transmission

other tools for reflections


Outlook
-------------------------------------

While this tool is built spesifically for DF-XRM, other techniques, such as Bragg scanning probe X-ray microscopy, Bragg coherent diffraction imaging and Bragg X-ray ptychography face similar challenges, and while the current implementation of `DF-XRM_vis` will be useful for these techniques as well, developing separate variants for the different tools may serve the community better. We also note that DF-XRM was first implemented at the European Synchrotron Radiation Facility in 20??, the Stanford Linear Accelrator in 2021, and the Advanced light Source in 202?. Implementations at other synchrotrons and X-ray free electron lasers are expected during this decade.


