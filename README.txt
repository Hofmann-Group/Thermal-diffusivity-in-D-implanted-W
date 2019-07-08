Hi ! Below is a brief description of what the folders contain and a guide of how the codes and data in them can be used to create the figures/data reported in the paper



Folders - 

EM images - contains the SEM images of key locations in all samples showing the blistering effects, also contains cross-sectional micrographs from the LTLD sample


TC SAW Profiles - this folder contains the raw TGS data, the processing/fitting code and the output data with thermal diffusivity for all samples. 'cleaned' profiles--
                - indicate that data points outside the region of interest are removed. These profiles are figure 3 in the main paper and the thermal diffusivity and SAW
                - profiles in the supplementary, given in the plots subfolder 


Fluence Temperature Profiles - this folder contains the raw temperature and fluence measurement data in a excel file. The dose_temp_plots_initial processes this data, --

                             -- applies BC's etc... and creates the data that is in the Processed Profiles folder. These processed profiles for all samples are then  -->
                             -- plot using the dose_temp_plots_paper to give figure 1 in the paper, given in the figures subfolder 

2D TC Temp Fluence Plot  - this folder contains the plotting_2D code and plotting_2D_2 codes. The plotting_2D code takes in the processed temperature and fluence profiles-
                         - and the processed TGS data, and generates the data for the variation of thermal diffusivity with temperature and fluence. plotting_2D dose this -
                         - for the high temperature exposure samples and plotting_2d_2 does it for the low temperature samples. The last two sections of plotting_2D then
                         - generates the final 2D plot that is figure 4 of the paper (2D_plot_TC_temp_fluence.fig)

Casino EDX Probing Depths - this folder contains the x-ray probing depths obtained for the EDX electron probing energies, using the CASINO Monte Carlo simulations
                          - the probing_depths code loads the depth_distributions of the x-rays and finds the 1/e depth, and plots the profiles too. 



See comments in codes for further information
contact: mohamed.reza@eng.ox.ac.uk or felix.hofmann@eng.ox.ac.uk for further info 
