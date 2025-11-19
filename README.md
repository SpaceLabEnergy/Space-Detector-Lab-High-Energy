# Space-Detector-Lab-High-Energy

### Purpose: 
Script to process data from three High-Energy Detectors and produce energy-channel calibration curves along with the absolute and intrinsic efficiencies and energy resolutions of each detector. 

### Important Operating Instructions:
Each detector has its own main .py script: **BGO.py**, **CdTe.py**, and **NaI.py**. When you run each script you will get outputs for that specific detector: if you run the BGO.py script you will only get results for the BGO detector. 

Additionally, important functions that these main codes call on are stored within the **gamma_ray_functions.py** script. The main scripts WILL NOT WORK if the **gamma_ray_functions.py** script is not also downloaded and kept within the same directory that the main script is run in. 

The main .py detector scripts (**BGO.py**, **CdTe.py**, and **NaI.py**) and the **gamma_ray_functions.py** script can be found in the **Scripts** folder

The data processed by the scripts can be found within the folder **Data for Code** and must be downloaded within the directory you run the main scripts in. 

To successfully run the code, the following libraries must be installed in the python environment:
- argparse
- numpy
- matplotlib.pyplot
- scipy
- lmfit

<ins>Suggested Procedure for Running the Code</ins>:
1. Create a folder/directory workspace
2. Download the scripts from the **Scripts** folder
3. Download the data from the **Data for Code** folder
4. Move the downloaded data and scripts into the created workspace (make sure there are no nested folders)
6. Verify your python environment has the required libraries downloaded
7. Run the code for each detector

<ins>Input</ins>: No initial user input, ensure the code is run within the same directory as the spectrum data within the **Data for Code** folder and the **gamma_ray_functions.py** script.

<ins>Output</ins>: Graphs of each element's spectrum with gaussian fits to photopeaks, graph of the detector's energy/channel calibration, graphs of the detector's intrinsic efficiency, absolute efficiency, and energy resolution as functions of peak energy, and a graph of Americium's efficiencies (intrinsic and absolute) as a function of off-angle response. 

