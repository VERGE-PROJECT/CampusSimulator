# Campus Scenario RAN Simulator
2025 - Mobile Communications Research Group - UPC

The software included here incorporates different Matlab scripts that simulate a Radio Access Network (RAN) in the Campus Nord scenario of Universitat Politècnica de Catalunya (UPC) where a number of base stations are deployed providing service to different User Equipments (UE). Some of these UEs can be configured to act as relays providing service to other UEs. 
The following figure illustrates the Campus Nord scenario. The environment is a 350 m x 125 m area with 25 buildings of 3 floors. The modelled area corresponds to the rectangle highlighted in red. The names of the buildings A1,...,D6 are also included in the figure. The positions of the three base stations are also included.
 
![image](https://github.com/user-attachments/assets/287a624c-8f7c-4b19-aac1-fb1a4597d8d7)

This campus scenario is the same where the measurements of the "Space/Time User Distribution Dataset" of the VERGE open dataspace have been collected (see https://github.com/VERGE-PROJECT/Space-Time-User-Distribution-Dataset). Then, the files included here could be modified to generate UE distributions in accordance with this dataset.

# DESCRIPTION OF THE MAIN SCRIPTS
The main scripts and functions included are:

•	campus_create_scenario.m: This script allows creating a scenario with a given number of buildings, streets, etc. Normally, this only needs to be executed once and it creates a file “CampusScenarioBuildings.mat” to be used in the different simulations. By default, this script is prepared to configure the scenario of the Campus Nord of UPC shown in previous figure.

•	campus_add_BSs.m: This script allows adding a number of BSs in the scenario and computing the propagation losses. This needs to be executed once for each configuration of BSs that we would like to test. It needs the CampusScenarioBuildings.mat to be previously created and it generates a file "CampusScenario_3BSs.mat”. By default, the script considers the three base stations illustrated in previous figure.

•	campus_add_relays.m: This script is used to generate the Relay Database. Each file in the relay database includes the propagation losses of a relay at a certain position in the scenario. It needs the CampusScenarioBuildings.mat to be previously created. The relay database only needs to be executed once for a given scenario and for a given configuration of the relay propagation model.

•	main.m: This is the script that executes the simulation runs specified by config.num_simulations, and generates the results file with the different statistics. Each simulation run corresponds to a different random distribution of the UEs acting as relays selected from the relay database. Different simulations for different numbers of relays can be executed by indicating the number of relays to simulate in the variable parameter_to_change.

•	campus_prop_model.m: computation of path loss for different types of propagation models

•	shadowing_2D.m: Generates a 2D map of correlated shadowing

•	mobility_model_pedestrian.m: Computes the position update for a pedestrian that moves through the streets of the scenario.

•	campus_mobility_regions.m: This script allows configuring the areas where the pedestrians can move.

The rest of included scripts are called by the ones listed above.

# EXECUTION PROCEDURE

Step 1) Create the scenario with the buildings and the trajectories:

- Configure the scenario size and the buildings positions at: campus_create_scenario.m
- Configure the streets for pedestrian and the pedestrian random walk areas at mobility_regions.m
- Execute script:  campus_create_scenario.m   (this already calls campus_mobility_regions script)
The resulting scenario is stored at:   CampusScenarioBuildings.mat

Step 2) Add the BSs and compute the propagation losses:
- Configure the BS positions and the parameters that impact the prop. model (freq., heights) at file campus_add_BSs.m
- Execute the script campus_add_BSs.m
The resulting scenario with the BSs is stored at: CampusScenario_3BSs.mat

Step 3) Create the relay database (this only requires step 1 to be executed before, but not necessarily step 2).
- Configure the relay parameters that impact on the propagation model at campus_add_relays.m
- Execute the script campus_add_relays.m
The resulting relays are stored in the specified directory (e.g.  ".\CampusRelayDatabase3.5GHz\")

Step 4) Simulation (this requires that the scenario file ‘CampusScenario_3BSs.mat’ is computed before).
- Configure the parameters of the simulation (e.g. Tx powers, noise, mobility model parameters, etc.)
- Execute script main.m 

# RESULTS
The execution of main.m generates the following files:

- Results file (name of the file is specified in "config.output_file_name" variable): Details of the statistics per simulation and in global. Also includes the configuration parameters. 
- map file (name specified in "config.map_BS_file_name" variable): Details of the BSs and maps of SNR_total, speff_total and serving BSs. This file is only generated if no file exists with the same name in the same directory. This allows that, if we are doing simulations just changing the number of UEs, mobility models, etc. this does not need to be resaved every time. If we are doing simulations changing the tx_power, antenna_gains, etc., it is better to give different names to the file to reflect the considered conditions.
