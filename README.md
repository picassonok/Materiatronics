# Materiatronics
Ready-to-use semi-automated analysis of arbitrary dipolar meta-atoms


This file contains information on how to use the modular decomposition introduced in original paper [V.S. Asadchy and S.A. Tretyakov, “Materiatronics: Modular analysis of arbitrary meta-atoms”, Physical Review Applied, 2019]. The folder, in which this file is located, for user’s convenience has already some example files corresponding to the decomposition of a split-ring resonator shown in Fig. 3a of the paper. Thus, following steps below, you can use the enclosed files or create files for your own meta-atom decomposition. See also the video file in the supplementary demonstrating all these steps.
1. Copy and paste your meta-atom into ANSYS HFSS (only in version 2017 and newer) project "Materiatronics_polarizabilities.aedt" (CTRL-C / CTRL-V will work). If you just want to use the enclosed files, just skip items 1.-6.
2.  Specify in the Design Properties the value for "min_freq" (the minimum frequency you will use for calculating the polarizabilities).
3.  Specify in "Analysis-Setup1" the frequency (usually, it is the maximum frequency used in your frequency sweep) and maximum number of adaptive passes.
4.  Specify the parameters in "Sweep" tab: "Start" and "End" frequencies, number of sweep points.
5.  Click button "Analyze All" to calculate the project.
6.  After calculations are ready, Export the results from plot "all_polarizabilities" (right click - Export) as a .csv file. Also check the units along the X and Y axes in this plot (e.g., [GHz] and [mV]), you will need them since ANSYS HFSS does not include this information into an exported file. 
7.  Open MATLAB file "Materiatronics.m". Update rows 6-9 if needed. Run the code.
8.  If your inclusion possesses strong spatial dispersion, the polarizability extraction technique used in this Matlab code may have some error, which in turn will lead to an error for the modular decomposition. The relative magnitude of this error is calculated in variable "relative_error". Typically, error below 0.3 (meaning 30%) provides acceptable results.
9.  The MATLAB code will automatically generate various polarizability plots as well as a new .py file with the same name as your initial .csv file. The code also creates a table with all the data of modular decomposition of the meta-atom. If you want to see a visual 3D illustration of the modular decomposition of your meta-atom, follow steps 10 and 11.
10. Run ANSYS HFSS file "Materiatronics_modules.aedt" (do not change its name). Follow commands "Tools - Run Script...", choose the newly created .py file and in about 50 sec enjoy the visual appearance of the decomposition. If you have not created your .py file yet, you can use instead the enclosed file ‘example_file.py’ to see how it works for the decomposition of a split-ring resonator. Note that ANSYS may crash once during the script running. In this case just run the project and script again. Modules with names “Max_chiral_orientation” and “Max_Tellegen_orientation” (in black) depict the orientations of maximum chiral and Tellegen bianisotropic couplings. 
11. In case you want to reset the decomposition in file "Materiatronics_modules.aedt" to the initial state, run the attached script with name "resetting the modules.py".

For any assistance with the computational code and modular decomposition, you can contact Dr. Viktar Asadchy (viktar.asadchy@gmail.com). 
