%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates a number of relays and stores it in the data base.
% It requires the following variables to be previously computed:
% - map_valid_points_relays: Pixels where a relay can be placed.
% - map_buildings: Map with the buildings of the scenario (each pixel
% indicates the building number)
% - map_indoor_points: Pixels indicating indoor points
% - buildings: Structure with the vertices of each building.
%
% All the variables of the scenario (sizeX,sizeY, pixel_size, etc.) should
% also have been previously defined.
%
% (c) 2025 - Mobile Communications Research Group - UPC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%% Requires the maps with the buildings to be previously generated.
clear all;
name_scenario_buildings_file='CampusScenarioBuildings.mat';
if exist(name_scenario_buildings_file)==0
    fprintf('\nSCENARIO BUILDINGS FILE NOT AVAILABLE. IT HAS TO BE GENERATED USING create_scenario.m');
    stop
else
    load(name_scenario_buildings_file);
end


%% 
%CONSTANTS:
D2D_Siemens=4;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%CONFIGURATION PARAMETERS (RELAY LINKS) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_relays=10; %%Number of relays to generate
directory_relay_database='CampusRelayDatabase3.5GHz\';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Propagation model of BS-UE link
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_D2D=D2D_Siemens;
f_relay=3.5;  %GHz
loss_per_floor=21; %dB (Considered when tx and rx of D2D link are in different floors)
sigma_LOS_relay=0;
sigma_NLOS_relay=12;
dcorr_LOS_relay=10; %The model does not specify any value.
dcorr_NLOS_relay=13; %The model does not specify any value.
sigma_penetration_relay=4.4; %dB Set the same value as for BSs
dcorr_penetration_relay=7; %m (3GPP TR 38.901 Table 7.5-6 - Correlation Distance SF)
resolution_dist_NLOS_relay=5; %Resolution when generating NLOS zones according to the NLOS probability (better put >1, otherwise the map will be the same for all relays)
dmax_NLOS_relay=50; %Beyond this distance it will be the NLOS.
dmax_computation_relay=200; %Maximum distance to compute the path loss of a relay
pixels_max=ceil(dmax_computation_relay/pixel_size);

PLmax=180; %dB Points with PL>PLmax are considered as 0 to save space



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTATIONS            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.- Generation of the relays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for n=1:num_relays
    if mod(n,100)==0
        fprintf('Generating Relay: %d\n',n);
    end
    %Choose a position
    rng(n);  %For reproducibility control the random position
    
    valid_position=0;
    %When selecting the position we need to make sure that it is placed in
    %the valid points (map_valid_points_relays)
    while valid_position==0
        %Repeat until finding a valid position of the relay
        relay.pos(1)=scenario_x*rand;
        relay.pos(2)=scenario_y*rand;
        posrelay_index=ceil(relay.pos/pixel_size);
        if map_valid_points_relays(posrelay_index(1),posrelay_index(2))==1
            valid_position=1;
        end
    end
    
    relay.building=map_buildings(posrelay_index(1),posrelay_index(2));
    if relay.building==0
        relay.floor=1;
    else
        %Indoor position: Choose randomly a floor:
        relay.floor=randi(num_floors);
    end
    relay.hUT_total=hUT+(relay.floor-1)*floor_height;
    
    map_NLOS_buildings = find_NLOS_points_Relay(sizeX,sizeY,pixel_size,relay.pos,buildings,map_buildings,dmax_NLOS_relay);
    map_NLOS_model = genera_map_NLOS(sizeX,sizeY,pixel_size,resolution_dist_NLOS_relay,dmax_NLOS_relay,relay.pos,model_D2D,hUT,f_relay,10*n);
    relay.map_NLOS=max(map_NLOS_buildings,map_NLOS_model);
    relay.map_indoor_distances=find_indoor_distances_relay(sizeX,sizeY,pixel_size,relay.pos,buildings,map_buildings,dmax_computation_relay);

    
    relay.map_PL{1}=1000*ones(sizeX,sizeY); %L->infinite
    for flo=2:num_floors
        relay.map_PL{flo}=1000*map_indoor_points; %L->infinite in all indoor points
    end
    
    for i=max(1,posrelay_index(1)-pixels_max):min(sizeX,posrelay_index(1)+pixels_max)
        for j=max(1,posrelay_index(2)-pixels_max):min(sizeY,posrelay_index(2)+pixels_max)
            pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
            dist=norm(pos-relay.pos); %2D distance
            if dist<=dmax_computation_relay    
                %%CONSIDER DIFFERENT POSSIBILITIES
                if relay.building==0 && map_buildings(i,j)==0
                    %Outdoor-to-outdoor:   (floor=1)
                    relay.map_PL{1}(i,j)=prop_model(dist,model_D2D,f_relay,relay.hUT_total,hUT,1,relay.map_NLOS(i,j),map_indoor_points(i,j),relay.map_indoor_distances(i,j));
                    %Note that when calling the prop_model function the hE
                    %input is not relevant for D2D model, so we set it at
                    %1.
                elseif (relay.building>0 && map_buildings(i,j)==0) 
                    %Relay indoor and point outdoor:
                    %Only compute for flo=1
                    dist3D=sqrt(dist^2+(relay.hUT_total-hUT)^2); %3D distance
                    %When computing map_PL we put the indoor condition to 1
                    %to include the outdoor-to-indoor losses.
                    relay.map_PL{1}(i,j)=prop_model(dist3D,model_D2D,f_relay,relay.hUT_total,hUT,1,relay.map_NLOS(i,j),1,relay.map_indoor_distances(i,j));
                    %Add the loss for the difference of floors.
                    relay.map_PL{1}(i,j)=relay.map_PL{1}(i,j)+loss_per_floor*(relay.floor-1);
                elseif (relay.building==0 && map_buildings(i,j)>0)
                    %Relay outdoor and point indoor.
                    %We compute for all the floors.
                    for flo=1:num_floors
                        dist3D=sqrt(dist^2+(relay.hUT_total-hUT_total(flo))^2); %3D distance
                        relay.map_PL{flo}(i,j)=prop_model(dist3D,model_D2D,f_relay,relay.hUT_total,hUT_total(flo),1,relay.map_NLOS(i,j),map_indoor_points(i,j),relay.map_indoor_distances(i,j));
                        %Add the loss for the difference of floors.
                        relay.map_PL{flo}(i,j)=relay.map_PL{flo}(i,j)+loss_per_floor*(flo-1);
                    end
                elseif relay.building==map_buildings(i,j)
                    %Relay indoor and point indoor inside the same building
                    %(different from 0)
                    %We compute for all the floors
                    for flo=1:num_floors
                        dist3D=sqrt(dist^2+(relay.hUT_total-hUT_total(flo))^2); %3D distance
                        %When computing the map_PL we put the indoor_distance
                        %equal to 0 and the indoor point as 0. In this way
                        %we don't take into account the indoor-to-outdoor
                        %losses.
                        relay.map_PL{flo}(i,j)=prop_model(dist3D,model_D2D,f_relay,relay.hUT_total,hUT_total(flo),1,relay.map_NLOS(i,j),0,0);
                        %Add the loss for the difference of floors.
                        relay.map_PL{flo}(i,j)=relay.map_PL{flo}(i,j)+loss_per_floor*abs(flo-relay.floor);
                    end
                end %If not condition is met, it means that the relay and the point are indoor but they are in different buildings, so the PL is just the default infinite value.
            end
        end
    end
    
    %%SHADOWING LOSSES:
    %Assumptions for shadowing:
    %Outdoor relay: same computations like for a BS
    %   -For indoor points, penetration loss is added on top of the shadowing
    %   -The shadowing loss is the same in all floors
    %   -The sigma of penetration losses is the same for relays as for BSs.
    
    relay.shadow_LOS=shadowing_2D(sizeX,sizeY,pixel_size,sigma_LOS_relay,dcorr_LOS_relay);
    relay.shadow_NLOS=shadowing_2D(sizeX,sizeY,pixel_size,sigma_NLOS_relay,dcorr_NLOS_relay);
    if relay.building==0    %Outdoor relay
        for flo=1:num_floors
            relay.shadow_penetration{flo}=shadowing_2D(sizeX,sizeY,pixel_size,sigma_penetration_relay,dcorr_penetration_relay);
            if flo==1
                relay.map_totalloss{flo}=relay.map_PL{flo}+relay.shadow_LOS.*(1-relay.map_NLOS)+relay.shadow_NLOS.*relay.map_NLOS+relay.shadow_penetration{flo}.*map_indoor_points;
            else
                %Only consider the indoor points (we multiply by the map of
                %indoor points)
                relay.map_totalloss{flo}=map_indoor_points.*(relay.map_PL{flo}+relay.shadow_LOS.*(1-relay.map_NLOS)+relay.shadow_NLOS.*relay.map_NLOS+relay.shadow_penetration{flo}.*map_indoor_points);
            end
        end
    else   %Indoor relay
        %Shadow penetration only counted for outdoor points in floor=1.
            
        for flo=1:num_floors    
            if flo==1
                relay.shadow_penetration{flo}=shadowing_2D(sizeX,sizeY,pixel_size,sigma_penetration_relay,dcorr_penetration_relay);
                %Shadow penetration only accounts for the outdoor points.
                relay.map_totalloss{flo}=relay.map_PL{flo}+relay.shadow_LOS.*(1-relay.map_NLOS)+relay.shadow_NLOS.*relay.map_NLOS+relay.shadow_penetration{flo}.*(1-map_indoor_points);
            else
                %No shadow penetration.
                %Total loss is only computed for the indoor points.
                relay.map_totalloss{flo}=map_indoor_points.*(relay.map_PL{flo}+relay.shadow_LOS.*(1-relay.map_NLOS)+relay.shadow_NLOS.*relay.map_NLOS);
            end    
        end

    end
    
    
    %RELAY STORAGE:
    relay_stored.floor=relay.floor;
    relay_stored.pos=relay.pos;
    for flo=1:num_floors
        relay_stored.map_totalloss{flo}=relay.map_totalloss{flo};
        relay_stored.map_totalloss{flo}(relay.map_totalloss{flo}>180)=0;
    end
    name_output_file=[directory_relay_database,'relay_',num2str(n),'.mat'];
    save(name_output_file,'relay_stored');   
end 

 
    
   
 


    
    
 
