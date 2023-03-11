Anyone is free to contribute to or fork DF-XRM_vis

**Would you like to add a new lens configuration or material to the dropdown menus?**
- Make a pull request, create an issue on github, or email tmara@dtu.dk

**Did you find a bug?**
- Please create an issue on github and/or email tmara@dtu.dk
- Alternatively: open a new GitHub pull request with a patch

**Do you wish to develop this code for a different experiment type?**
- I would advise developing a new geometry in a forked repository. 
- Email tmara@dtu.dk if you want an orientation of the current codebase.

**Do you want to utilize part of this code in a bigger project (for example together with simulations)?**
- You are free to use any part of of this codebase for any purpose as outlined in the lisence. 
- Please cite the DF-XRM_vis if you found it useful.

**Other ideas for how to improve DF-XRM_vis?**
- Create an issue on github, a pull request or email tmara@dtu.dk

*On automated tests:*
- DF-XRM_vis employs an automated test (tests/test_geom.py) to prevent unitended changes to the geometry. 
- There are numerical checks to the angles and Q vector.
- There are also comparisons to template images. 
- Finally a pdf is generated.
- Please consult the pdf and look for changes when you push.
- If you make any changes the geometry or rendering of the figures, make sure you upload new versions of 2dfig.png, 3dfig.png and/or fig_optics.png in tests/, otherwise the test will fail.
- The figures are generated (in the root directory) by running pytest locally.
- Certain changes to DF-XRM_vis_streamlit.py will need to be reflected in tests/test_geom.py
