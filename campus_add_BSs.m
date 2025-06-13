%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script adds the BSs and generates the maps of total losses per BS
%
% (c) 2025 - Mobile Communications Research Group - UPC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Requires the maps with the buildings to be previously generated.
clear all;
rng(700);
name_scenario_buildings_file='CampusScenarioBuildings.mat';
if exist(name_scenario_buildings_file)==0
    fprintf('\nSCENARIO BUILDINGS FILE NOT AVAILABLE. IT HAS TO BE GENERATED USING create_scenario.m');
    stop
else
    load(name_scenario_buildings_file);
end

%% 
%CONSTANTS:
UMi=0;
UMa=1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Parameters of the propagation model of BS-UE link
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model=UMa;
hE=1;   %m   (valid only for UMi)
sigma_LOS=4;
dcorr_LOS=37;
sigma_NLOS=6;
dcorr_NLOS=50;
sigma_penetration=4.4; %dB (Low Loss model: 4.4 dB, High Loss model: 6.5 dB)
dcorr_penetration=7; %m (3GPP TR 38.901 Table 7.5-6 - Correlation Distance SF)
resolution_dist_NLOS=20;  %Resolution when generating NLOS zones according to the NLOS probability
dist_max_NLOS=sqrt(scenario_x^2+scenario_y^2); %Beyond this distance the model will be NLOS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of the BS (transmitter) 
%(ONLY THOSE THAT ARE RELEVANT FOR THE PROPAGATION MODEL. 
% THE OTHERS (POWER, BW, ETC. ARE SET IN THE SIMULATION SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=3.72576; %GHz
hBS = zeros(3,1);
hBS(1)=27; %m
hBS(2)=8; %m
hBS(3)=23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the UE (receiver)
%(ONLY THOSE THAT ARE RELEVANT FOR THE PROPAGATION MODEL. 
% THE OTHERS (POWER, BW, ETC. ARE SET IN THE SIMULATION SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No parameters.
% Note that the receiver height is specified in the "create_scenario.m".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Positions of the BSs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BS(1).pos=[-88,-43];  %m % Bruc
BS(2).pos=[243,182];     % Omega
BS(3).pos=[360,-170];    % Esade

numBS=size(BS,2);

intersite_distances=zeros(numBS,numBS);
min_ISD=1E6+zeros(numBS,1);
for n=1:numBS
    for m=1:numBS
        if n~=m
            intersite_distances(n,m)=norm(BS(m).pos-BS(n).pos);
            if intersite_distances(n,m)<min_ISD(n)
                min_ISD(n)=intersite_distances(n,m);
            end
        end 
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.- Computation of NLOS zones for each BS according to the 
%buildings and the LOS probabiliy model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_aleat_BS(1)=8;
num_aleat_BS(2)=10;
num_aleat_BS(3)=22;
for n=1:numBS
    fprintf('Computing NLOS map for BS:%d \n',n);
    map_NLOS_model = genera_map_NLOS(sizeX,sizeY,pixel_size,resolution_dist_NLOS,dist_max_NLOS,BS(n).pos,model,hUT,f,num_aleat_BS(n));
    map_NLOS_buildings = find_NLOS_points(sizeX,sizeY,pixel_size,BS(n).pos,buildings,map_indoor_points);
    BS(n).map_NLOS=max(map_NLOS_model,map_NLOS_buildings);
    BS(n).map_indoor_distances=find_indoor_distances(sizeX,sizeY,pixel_size,BS(n).pos,buildings,map_buildings);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.- Computation of propagation losses per BS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:numBS
    fprintf('Computing Pathloss for BS:%d \n',n);
    for flo=1:num_floors
        BS(n).map_PL{flo}=zeros(sizeX,sizeY);
    end
    for i=1:sizeX
        for j=1:sizeY
            pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
            dist=norm(pos-BS(n).pos);
            for flo=1:num_floors
                BS(n).map_PL{flo}(i,j)= campus_prop_model(dist,model,f,hUT_total(flo),hBS(n),hE,BS(n).map_NLOS(i,j),map_indoor_points(i,j),BS(n).map_indoor_distances(i,j));     
                if map_indoor_points(i,j)==0
                    %Outdoor point. Only compute for flo=1, so we can break
                    %the loop
                    break
                end
            end
        end
    end
    

    %Assumptions for shadowing:
    %-For indoor points, penetration loss is added on top of the shadowing
    %-The shadowing loss is the same in all floors
    
    BS(n).shadow_LOS = shadowing_2D(sizeX,sizeY,pixel_size,sigma_LOS,dcorr_LOS);
    BS(n).shadow_NLOS = shadowing_2D(sizeX,sizeY,pixel_size,sigma_NLOS,dcorr_NLOS);
    
    

   
    
    for flo=1:num_floors
       BS(n).shadow_penetration{flo}=shadowing_2D(sizeX,sizeY,pixel_size,sigma_penetration,dcorr_penetration);
       if flo==1
           BS(n).map_totalloss{flo}=BS(n).map_PL{flo}+BS(n).shadow_LOS.*(1-BS(n).map_NLOS)+BS(n).shadow_NLOS.*BS(n).map_NLOS+BS(n).shadow_penetration{flo}.*map_indoor_points;
       else
           %Only consider the indoor points (we multiply by the map of
           %indoor points)
           BS(n).map_totalloss{flo}=map_indoor_points.*(BS(n).map_PL{flo}+BS(n).shadow_LOS.*(1-BS(n).map_NLOS)+BS(n).shadow_NLOS.*BS(n).map_NLOS+BS(n).shadow_penetration{flo}.*map_indoor_points);
       end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save scenario with the BS maps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%To save space we keep only the relevant maps in the BS struct
for n=1:numBS
   BS_to_store(n).pos=BS(n).pos;
   BS_to_store(n).map_totalloss=BS(n).map_totalloss;
end
BS=BS_to_store;
clear BS_to_store;

save('CampusScenario_3BSs.mat','-v7.3');
