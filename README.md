# BEAMS
 BEAMS processing code

Code to extract data from the SeapHOx with two pumps in the BEAMS configuration. Deployed on a reef, the BEAMS pumps sea water from two heights to measure the carbon flux over a reef. Combined with an accoustic doppler and a PAR sensor, this information can be used to estimate productivity (NEP and NEC). 

**DEPENDENCIES**
Clone the following repositories and have in MATLAB path: 
https://github.com/CarbonLab/import_functions/tree/main/lab_instruments/pHETSOS
https://github.com/yui-takeshita/seawater_equilib

Download the following functions and have in MATLAB path:
https://www.marine.csiro.au/datacentre/projects/seawater/seawater_ver1_2d.zip
https://www.mathworks.com/matlabcentral/fileexchange/23441-polar-cartesian-into-degree-north-reference


**HOW TO RUN THE SCRIPT**
"Extract_SP_BEAMS.m" is the wrapper script to extract data and calculate all other parameters.
Inputs are defined in the first section of the script under Meta information.

base_dir: base directory where the BEAMS (SeapHOx) and aquadopp data are stored. Edit to match the location on local machine.

mat_name: file save location. Set up to save in the base directory, can be changed to the save in a local folder.

use_def_Sal: Can set to true if the MCAT is broken and have to use an estimated value form bottle samples
startdate/enddate: Manually added indices to trim the data to the desired start and end times. MUST FILL OUT PROPERLY BEFORE RUNNING SCRIPT.

TA0: Background total alkalinity can be adjusted/tuned for location of deployment based on bottle samples.
Q: photosynthetic quotient. Can also be tuned for deployment location.

Sal_est: Simple estimate to give a model a starting location.

hasAquadopp: Default is true.

delimit: Delimiter of the csv output from the Aquadopp. In most cases it is ';'
min2avg: can be tweaked but is generally set to the BEAMS sampling frequency.

pumpz: Need more info on this one but should leave be for now.
ipumps: Pump IDs to use for calculations. Most BEAMS set ups only use pump options 1 and 2.

meta_parser: function to parse through the meta output from GetMfetData. Not a bullet proof solution and will likely be updated in the future.

