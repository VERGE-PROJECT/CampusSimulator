%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates the scenario (buildings, streets, trajectories)
%
% (c) 2025 - Mobile Communications Research Group - UPC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
rng(700);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONFIGURATION PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Scenario parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenario_x=335;
scenario_y=125;
pixel_size=1;
sizeX=scenario_x/pixel_size;
sizeY=scenario_y/pixel_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Receiver height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hUT=1.5;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Characterisation of the buildings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_floors=3; %Note: floor=1 is ground.
floor_height=3.5; %m  Note that the models are valid for hUT<=22.5m. We should take
%this into account when setting floor_height and number of floors.
for flo=1:num_floors  %Note: we cannot use the name "floor" for a variable because there is a Matlab function called floor.
    hUT_total(flo)=hUT+(flo-1)*floor_height;
end

num_buildings=24; 
% -- D1-D6 --
buildings(1).vertices=[5,6; 107,6; 107,25; 5,25]; %Coordinates in m
buildings(2).vertices=[112,6; 151,6; 151,25; 112,25];
buildings(3).vertices=[156,6; 195,6; 195,25; 156,25];
buildings(4).vertices=[200,6; 239,6; 239,25; 200,25];
buildings(5).vertices=[244,6; 283,6; 283,25; 244,25];
buildings(6).vertices=[288,6; 327,6; 327,25; 288,25];
% -- C1-C6 --
buildings(7).vertices=[18,35; 107,35; 107,55; 18,55; 18,41];
buildings(8).vertices=[112,35; 151,35; 151,55; 112,55];
buildings(9).vertices=[156,35; 195,35; 195,55; 156,55];
buildings(10).vertices=[200,35; 239,35; 239,55; 200,55];
buildings(11).vertices=[244,35; 283,35; 283,55; 244,55];
buildings(12).vertices=[288,35; 327,35; 327,55; 288,55];
% -- B1-B6 --
buildings(13).vertices=[50,67; 90,67; 90,86; 50,86];
buildings(14).vertices=[129,67; 151,67; 151,86; 129,86];
buildings(15).vertices=[156,67; 178,67; 178,86; 156,86];
buildings(16).vertices=[217,67; 239,67; 239,86; 217,86];
buildings(17).vertices=[244,67; 266,67; 266,86; 244,86];
buildings(18).vertices=[306,67; 327,67; 327,86; 306,86];
% -- A1-A6 --
buildings(19).vertices=[68,100; 107,100; 107,119; 68,119];
buildings(20).vertices=[112,100; 151,100; 151,119; 112,119];
buildings(21).vertices=[156,100; 195,100; 195,119; 156,119];
buildings(22).vertices=[200,100; 239,100; 239,119; 200,119];
buildings(23).vertices=[244,100; 283,100; 283,119; 244,119];
buildings(24).vertices=[288,100; 327,100; 327,119; 288,119];


num_buildings=size(buildings,2);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTATIONS            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.- Computation of the positions that are indoor and the map 
% of buildings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Computing Indoor map.\n');
[map_indoor_points,map_buildings]=find_indoor_points(sizeX,sizeY,pixel_size,buildings);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.- Create the mobility regions for pedestrians and cars
%and the valid points where to place relays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Creating Mobility Regions.\n');
campus_mobility_regions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.- Save scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('CampusScenarioBuildings.mat');