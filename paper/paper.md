---

author:
- 'Trygve Magnus Ræder, tmara\@dtu.dk'
date: April 2022
title: 'DF-XRM\_viz'
---



# Summary

`DF-XRM_vis` is a visualization tool to aid in planning, execution, and data analysis of X-ray experiments based on bragg diffraction at synchrotrons and X-ray free electron lasers.
The toolkit is built in python and has a streamlit interface, as well as a jupyterlab example.
The highlight of the toolkit is generating three dimensional models of experimental geometry, rendered together with the unit cell of the material. By hosting the streamlit application on a publically avaiable adress, the toolkit becomes available also to researchers without programming experience.

# Statement of need

Advances in X-ray brilliance and optics has recently led to the development of Bragg scanning probe X-ray microscopy, Bragg coherent diffraction imaging, Bragg X-ray ptychography and dark-field X-ray microscopy. 
All the aforementioned techniques fundamentally capture the same crystallographic information, and have different advantages and disadvantages.
As all the aforementioned techniques are built upom Bragg diffraction, they all require the instrument to be aligned so that Braggs law is fulfilled. 
Failure to apprehend or communicate the constraints has been a liability in the practical application of DF-XRM.
`DF-XRM_vis` is designed to be an aid in assessing the geometrical constraints that determine Bragg diffraction. `DF-XRM_vis` may also aid in data analysis, serve as a
teaching tool, or an aid when discussing the practical implications of advanced X-ray techniqes more broadly.

A successful experiment requires that the following pitfalls be avoided:

1. The sample absorption at the desired X-ray energy is too high and/or the sample is too thick.

2. The desired reflection and X-ray energy does not produce a diffraction angle suitable for the instrument.

3. The sample mount assumes a different crystallographic orientation than what the crystal actually has.

Avoiding the above pitfalls requires good communication of 3d geometry
between the researchers involved in the experiment. Traditionally, this
has been an (often challenging) exercise involving the use of a
whiteboard, but as more communication has shifted online, collaboration
on 3d geometry has become more difficult. I have personally observed
experiments loose days of beamtime to the issues listed above, as a result
of poor online communication in the lead-up to experiments.


While `DF-XRM_vis` is designed to address all the above points, we note that the novelty of `DF-XRM_vis` lies with adressing the 3rd point, which has not previously been adressed in any comparable tool. 
We will here describe the capabilities of the streamlit interface to `DF-XRM_vis` in order of increasing user input, which is the same as the order as the points listed above.

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

This tool is built initially for DF-XRM, and therefore inclueds the option of visualizing the lens spesific to this technique. Spesific options relating to other techniqes may be implemented in the future. We also note that DF-XRM was first implemented at the European Synchrotron Radiation Facility in 20??, the Stanford Linear Accelrator in 2021, and the Advanced light Source in 202?. Implementations at other synchrotrons and X-ray free electron lasers are expected during this decade. 


