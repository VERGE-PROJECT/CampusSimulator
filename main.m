%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Script that executes the SIMULATION PROCESS     %%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) 2025 - Mobile Communications Research Group - UPC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for parameter_to_change=[2,3] %Number of relays 

% parameter_to_change --> num_relays

clearvars -except parameter_to_change

name_scenario_file='CampusScenario_3BSs.mat';
if exist(name_scenario_file)==0
    fprintf('\nSCENARIO FILE NOT AVAILABLE. IT HAS TO BE GENERATED');
    stop
else
    load(name_scenario_file);
end
activate_debug_UE=0;
testUE=100; %Id of the UE for debugging purposes.
test_s=2; %Id of the group of testUE.

config.store_samples_in_file=0;  %If set to 0 it will remove the samples of spectral efficiency, connection times, etc. before storing the file (to prevent a very large file).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%CONFIGURATION PARAMETERS (For the simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config.min_speff_relay=1; %Minimum SpEff between the BS and the relay to allow being a relay.
config.allow_multihop=0;    %(0: relays only connected to the BS. 1: relays can be connected to other relays.
config.directory_relay_database='CampusRelayDatabase3.5GHz\';

config.num_simulations=5;  %Number of simulations to be executed. Each simulation corresponds to a different distribution of relays.
config.simulation_duration=1000;   %s
config.time_step=1; %s Duration of the simulation time step
config.output_file_name=['res_campus_indoor_P22_Smin1_R100M_f35_N',num2str(parameter_to_change),'_relays.mat'];
config.map_BS_file_name='CampusMap3BSs_P22.mat';  %Change the name of the file if executing simulations by varying the txpower parameters.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of the BS (transmitter) 
% (The parameters that affect the prop. model (like freq. and heights) are set in the
% create_scenario script).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config.gain_BS=26; %26  dB  (antenna BS)
config.tx_power(1)=22;  %25    dBm (ref: Table 7.8-1 of TR38.901)  (old value: 21 dBm)
config.tx_power(2)=22;
config.tx_power(3)=22;
config.bandwidth=100;  %MHz
config.tx_power_per_Hz(1)=config.tx_power(1)-10*log10(config.bandwidth*1E6); %dBm/Hz
config.tx_power_per_Hz(2)=config.tx_power(2)-10*log10(config.bandwidth*1E6); %dBm/Hz
config.tx_power_per_Hz(3)=config.tx_power(3)-10*log10(config.bandwidth*1E6); %dBm/Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of the relay (transmitter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config.gain_antenna_relay= 0;%3.0; %dB
config.tx_power_relay=-10;  %dBm
config.bandwidth_relay=100;  %MHz
config.tx_power_per_Hz_relay=config.tx_power_relay-10*log10(config.bandwidth_relay*1E6); %dBm/Hz
%Note that the frequency of the relay is defined with the propagation model parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of the UE (receiver)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For the direct link with the BS:
config.gain_UE=10; %dB (antenna UE)
config.noise_figure=9; %dB
config.noise_density=-174+config.noise_figure;  %dBm/Hz
config.speff_max=7.4063; %Maximum spectral efficiency (b/s/Hz)
config.speff_min=0.2344; %Minimum spectral efficiency (below this value it will be 0)

%For the relay link:
config.gain_antenna_UE_relay= 0;%3.0; %dB

config.speff_limit_outage=1; %b/s/Hz  - Min. speff for outage probability computations. (it could be equal to the variable "speff_min" or to another value.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of the mobility model (pedestrian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config.avg_time_direction=10; %s  %Average time that a pedestrian follows the same direction with Random Walk
config.prob_change_RW=1/(config.avg_time_direction/config.time_step); %Probability of changing direction in a position update
config.angle_change_RW=45*pi/180; %rad    %Range of angle when changing direction with RW (+angle_change_RW/-angle_change_RW)
config.prob_change_inters_ped=0.05;        %Probability of changing direction in an intersection point of a street.
config.speed_pedestrian=3/3.6;    %m/s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numbers of users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config.num_users=100;

config.num_stationary=80;
config.num_pedestrians=20;

% config.relay_involvement_factor=parameter_to_change/100; %Fraction of stationary UEs that act as relays

config.max_relays_database=size(dir(config.directory_relay_database),1)-2;  %Maximum number of available relays in the database (it's the number of files in the directory minus 2 that correspond to the . and ..)
% config.num_relays=min(ceil(config.num_users*config.ratio_stationary*config.relay_involvement_factor),config.max_relays_database);
config.num_relays=parameter_to_change;

config.session_rate_per_UE=1;  %sessions/h
config.avg_session_duration=300;       %s

config.lambda=zeros(1,3); %s=1:Pedestrians, s=2:stationary non relay, s=3: stationary relay, 

config.lambda(1)=config.num_pedestrians*config.session_rate_per_UE/3600;   %sessions/s
config.lambda(2)=config.num_stationary*config.session_rate_per_UE/3600;   %sessions/s
config.lambda(3)=1E-12;   %sessions/s;   %sessions/s
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.- Computation of SNR and SpEff maps per BS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:numBS
    fprintf('Computing SNR and SpEff map for BS:%d \n',n);
    %Computation of SNR and spectral efficiency
    for flo=1:num_floors
        BS(n).map_SNR{flo}=config.tx_power_per_Hz(n)+config.gain_BS+config.gain_UE-BS(n).map_totalloss{flo}-config.noise_density;
        BS(n).map_speff{flo}=min(log2(1+power(10,0.1*BS(n).map_SNR{flo})),config.speff_max);
        BS(n).map_speff{flo}(BS(n).map_speff{flo}<config.speff_min)=0;
        if flo>1
            BS(n).map_speff{flo}=map_indoor_points.*BS(n).map_speff{flo};
            BS(n).map_SNR{flo}=map_indoor_points.*BS(n).map_SNR{flo};
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.- Computation of total maps for the direct link
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Computing total maps.\n');
for flo=1:num_floors
    SNR_total_BS{flo}=zeros(sizeX,sizeY)-1000;
    speff_total_BS{flo}=zeros(sizeX,sizeY);
    serving_BS{flo}=zeros(sizeX,sizeY);
end

for n=1:numBS
    for flo=1:num_floors
        mask_BS_service_area=zeros(sizeX,sizeY);
        if flo==1
           mask_BS_service_area(BS(n).map_SNR{flo}>SNR_total_BS{flo})=1; %Map with the points where the UE will connect to the BS. We use the SNR instead of sp_eff because the max. SNR is not bounded.
        else
           mask_BS_service_area(BS(n).map_SNR{flo}>SNR_total_BS{flo})=1; %Map with the points where the UE will connect to the BS. We use the SNR instead of sp_eff because the max. SNR is not bounded. 
           mask_BS_service_area=mask_BS_service_area.*map_indoor_points; %Only consider the indoor points.
        end
        
        %Update of the total maps:
        SNR_total_BS{flo}=SNR_total_BS{flo}.*(1-mask_BS_service_area)+BS(n).map_SNR{flo}.*mask_BS_service_area;
        speff_total_BS{flo}=speff_total_BS{flo}.*(1-mask_BS_service_area)+BS(n).map_speff{flo}.*mask_BS_service_area;
        serving_BS{flo}=serving_BS{flo}.*(1-mask_BS_service_area)+n*mask_BS_service_area;
    end    
end   


%%
for sim=1:config.num_simulations
    fprintf('\nStarting simulation:%d',sim);
    rng(sim);
    
    for s = 1:3   
        for nBS = 1:numBS
            this_BS = ['BS_', num2str(nBS)];
            UEs_BS.(this_BS){s}.direct.indoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_BS.(this_BS){s}.direct.outdoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_BS.(this_BS){s}.indirect.indoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_BS.(this_BS){s}.indirect.outdoor = zeros(1,config.simulation_duration/config.time_step);
        end
        for nRelay = 1:config.num_relays
            UEs_relays{nRelay}.users{s}.indoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_relays{nRelay}.users{s}.outdoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_relays{nRelay}.users_outage{s}.indoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_relays{nRelay}.users_outage{s}.outdoor = zeros(1,config.simulation_duration/config.time_step);
            UEs_relays{nRelay}.avg_speff_improvement = 0; 
            UEs_relays{nRelay}.n_appearances = 0;
            UEs_relays{nRelay}.fraction_outage_users{s} = 0;
            
        end
    end
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 3) READ THE RELAYS MAPS AND COMPUTE THE SNR AND SPEFF
% % Resulting maps:
% %  -SNR_total_relay{flo} : Map of the SNR of the best link (relay or
% %  direct) in each point and floor. When connecting through the relay, the SNR is
% %  the minimum between the SNR in BS-Relay and Relay-UE
% %  -speff_total_relay{flo}: Map of the speff of best link (relay or direct)
% %  in each point/floor
% %  - serving_BS_relay_total{flo}: Map of the identifier of the BS or relay
% %  serving each point/floor. In case of relay n it is encoded as 10000+n.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %First we need to determine where is each relay connected (i.e. to a BS or
% %to another relay). In this way we can obtain the SNR and SPeff of the
% %BS-relay links and use this to compute the total maps.
% 
list_selected_relays=zeros(1,config.max_relays_database);
% %Select num_relays randomly from the number of relays in the database
% 
n=0; %Number of valid relays
while n<config.num_relays
    num_pending_relays=config.num_relays-n;
    list_candidate_relays=find(list_selected_relays==0);
    relay_numbers=randperm(size(list_candidate_relays,2),num_pending_relays);    

    %Read first the files to get the positions of all relays.
    for naux=1:size(relay_numbers,2)
        name_input_file=[config.directory_relay_database,'relay_',num2str(relay_numbers(naux)),'.mat'];
        load(name_input_file);
        pos_relay_index=[ceil(relay_stored.pos(1)/pixel_size),ceil(relay_stored.pos(2)/pixel_size)];
        
        if speff_total_BS{relay_stored.floor}(pos_relay_index(1),pos_relay_index(2))>=config.min_speff_relay
            %Valid relay:
            list_selected_relays(relay_numbers(naux))=1;
            n=n+1;
            relay_list(n).file_name=name_input_file;
            
            UEs_relays{n}.id = relay_numbers(naux);
            relay_list(n).floor=relay_stored.floor;
            %fprintf('Relay:%d, Floor:%d\n',n,relay_stored.floor);
            relay_list(n).pos=relay_stored.pos;
            relay_list(n).pos_relay_index=pos_relay_index;
    
            relay_list(n).serving_BS=serving_BS{relay_list(n).floor}(pos_relay_index(1),pos_relay_index(2));
            relay_list(n).speff_BS_relay=speff_total_BS{relay_list(n).floor}(pos_relay_index(1),pos_relay_index(2)); %SPeff of link BS-Relay (or link Relay-Relay)
            relay_list(n).SNR_BS_relay=SNR_total_BS{relay_list(n).floor}(pos_relay_index(1),pos_relay_index(2));   %SNR of link BS_relay (or link Relay-Relay)
            relay_list(n).serving_relay=0;   %This will be >0 in case the relay is served through another relay.
            relay_list(n).SNR_direct_link=relay_list(n).SNR_BS_relay;
            relay_list(n).speff_direct_link=speff_total_BS{relay_list(n).floor}(pos_relay_index(1),pos_relay_index(2));
        end
    end
    if n<config.num_relays && sum(list_selected_relays)==config.max_relays_database
       fprintf('\nTHERE ARE NOT SUFFICIENT VALID RELAYS, but process continues');
       config.num_relays=n;
    end
end

fprintf('\nRelay files read. Num_relays:%d\n',config.num_relays);
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%4.- OBTENTION OF THE TOTAL SNR AND SPEFF MAPS INCLUDING THE RELAYS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR_total_relay=SNR_total_BS;  %Map of SNR of the best link (relay or direct) in each point/floor. Initialised to the SNR of BSs.
speff_total_relay=speff_total_BS; %Map of speff of best link (relay or direct) in each point/floor. Initialised to the SNR of BSs.
serving_BS_relay_total=serving_BS; %Map of the identifier of the BS or relay serving each point/floor.

%num_map=1;
%stats.map_speff_total{num_map}=speff_total_relay;
%stats.map_serving_BS_total{num_map}=serving_BS_relay_total;

%%READ THE RELAY MAPS AND COMPUTATION OF SNR AND SPEFF
%We keep structure relay_list(n) where we store the position and floor only of
%the relays (not the maps, to save memory).
%num_relays=100;

for n=1:config.num_relays
    load(relay_list(n).file_name);
    if mod(n,100)==0
        fprintf('\nComputing map for relay: %d',n);
    end
    for flo=1:num_floors
        mask=zeros(sizeX,sizeY);
        mask(relay_stored.map_totalloss{flo}>0)=1; %These are the only valid points of the relay.
        
        aux.map_SNR{flo}=mask.*(config.tx_power_per_Hz_relay+config.gain_antenna_relay+config.gain_antenna_UE_relay-relay_stored.map_totalloss{flo}-config.noise_density);
        map_speff_ini=mask.*log2(1+power(10,0.1*aux.map_SNR{flo}));

        %The max. speff of the relay-UE link is bounded by the speff in the
        %BS-relay link.
        aux.map_speff{flo}=min(map_speff_ini,relay_list(n).speff_BS_relay);
        aux.map_speff{flo}(aux.map_speff{flo}<config.speff_min)=0;
        
        aux.map_SNR{flo}=min(aux.map_SNR{flo},relay_list(n).SNR_BS_relay); %We also bound the SNR.
        
        mask_relay_service_area=zeros(sizeX,sizeY);
        mask_relay_service_area(aux.map_SNR{flo}>SNR_total_relay{flo})=1; %Map with the points where the UE will connect to the relay. We use the SNR instead of sp_eff because the max. SNR is not bounded.
        mask_relay_service_area=mask.*mask_relay_service_area; %We need to multiply by the mask of valid points to avoid selecting invalid points.
        
        %Update of the total maps:
        SNR_total_relay{flo}=SNR_total_relay{flo}.*(1-mask_relay_service_area)+aux.map_SNR{flo}.*mask_relay_service_area;
        speff_total_relay{flo}=speff_total_relay{flo}.*(1-mask_relay_service_area)+aux.map_speff{flo}.*mask_relay_service_area;
        serving_BS_relay_total{flo}=serving_BS_relay_total{flo}.*(1-mask_relay_service_area)+(10000+n)*mask_relay_service_area;
        %Note: The points served by the relay are indicated with the number
        %10000+n.
    end
    %if mod(n,10)==0
        %Store map
    %    num_map=num_map+1;
    %    stats.map_speff_total{num_map}=speff_total_relay;
    %    stats.map_serving_BS_total{num_map}=serving_BS_relay_total;
    %end
end
fprintf('Reading relay maps completed.\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%5) DYNAMIC SIMULATION PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
for s=1:3
    num_active_UEs(s,1)=0;
    time_next_session_arrival(s,1)=(-1/config.lambda(s))*log(1-rand);
    num_generated_UEs(s,1)=0;
end
list_active_relays=zeros(config.num_relays,1);

stats(sim).num_active_UEs=zeros(ceil(config.simulation_duration/config.time_step),4);
stats(sim).avg_session_duration=zeros(1,4);
stats(sim).num_sessions_released=zeros(1,4);
stats(sim).num_samples_speff=zeros(1,4);
max_num_samples_estimation=ceil(max(config.lambda*config.avg_session_duration*2*ceil(config.simulation_duration/config.time_step))); %Preallocate a worst case size assuming in every time step 2*avg number of users.
stats(sim).samples_speff=zeros(max_num_samples_estimation,4);

stats(sim).num_samples_speff_indoor=zeros(1,2);  %Note: will be stored only for stationary users (s=2 or s=3)
stats(sim).samples_speff_indoor=zeros(max_num_samples_estimation,2);
stats(sim).num_samples_speff_outdoor=zeros(1,2); %Note: will be stored only for stationary users (s=2 or s=3)
stats(sim).samples_speff_outdoor=zeros(max_num_samples_estimation,2);


max_num_samples_connection_time=ceil(max(config.lambda*config.simulation_duration*20)); %Preallocate a worst case size assuming 20 changes of relay/BS per session.
stats(sim).num_samples_relay_connection_time=zeros(1,4);
stats(sim).samples_relay_connection_time=zeros(max_num_samples_connection_time,4);
stats(sim).num_samples_BS_connection_time=zeros(1,4);
stats(sim).samples_BS_connection_time=zeros(max_num_samples_connection_time,4);

stats(sim).num_samples_connect_relay=zeros(1,4);
stats(sim).num_samples_connect_relay_indoor=zeros(1,2);
stats(sim).num_samples_connect_relay_outdoor=zeros(1,2);



if config.num_relays>0
   stats(sim).relay_list=relay_list;  %Store the list of relays for each simulation.
end

time_index=0;
for time=0:config.time_step:config.simulation_duration
    time_index=time_index+1;
    if activate_debug_UE==0 && mod(time,100)==0
        fprintf('Sim:%d, Simulating time: %f Ped:%d, Stat:%d, Relay:%d\n',sim,time,num_active_UEs(1),num_active_UEs(2),num_active_UEs(3));
    end
    
    %MOBILITY OF ACTIVE UEs OF TYPE PEDESTRIAN 
    s=1; %Pedestrian UEs
    for i=1:num_active_UEs(s)
        [UElist{s}{i}.pos,UElist{s}{i}.direction,UElist{s}{i}.trajectory] = mobility_model_pedestrian(UElist{s}{i}.pos,UElist{s}{i}.direction,UElist{s}{i}.trajectory,UElist{s}{i}.speed,config.time_step,config.prob_change_RW,config.angle_change_RW,config.prob_change_inters_ped,sizeX,sizeY,pixel_size,map_valid_points_pedestrian,map_pedestrian_points,trajectories_ped,map_numbers_traject_ped);
    end
    
    
    %CHECK SESSION FINALISATIONS:
    for s=1:3
        if num_active_UEs(s)>0
            end_process=0;
            i=1;
        else
            end_process=1;
        end
        while ~end_process
            if UElist{s}{i}.end_session_time<=time
                %Before removing the UE, if it is a relay UE, let's
                %deactivate it from the list of active relays.
                if s==3
                    list_active_relays(UElist{s}{i}.associated_relay)=0;
                end
                stats(sim).num_sessions_released(s)=stats(sim).num_sessions_released(s)+1;
                stats(sim).avg_session_duration(s)=stats(sim).avg_session_duration(s)+UElist{s}{i}.end_session_time-UElist{s}{i}.start_time;
                
                %Measure statistic of connection time in BS or relay:
                if UElist{s}{i}.serving_relay>0
                   %It was connected to a relay:
                   stats(sim).num_samples_relay_connection_time(s)=stats(sim).num_samples_relay_connection_time(s)+1;
                   stats(sim).samples_relay_connection_time(stats(sim).num_samples_relay_connection_time(s),s)=UElist{s}{i}.end_session_time-UElist{s}{i}.start_time_serving_BS_or_relay;
                   if activate_debug_UE==1 && s==test_s && UElist{s}{i}.UEid==testUE
                        fprintf('\nTime:%d, UEid:%d STOPS. SampleTimeRelay:%4.2f',time,testUE,stats(sim).samples_relay_connection_time(stats(sim).num_samples_relay_connection_time(s),s));
                   end
                else
                   %It was connected to a BS:
                   stats(sim).num_samples_BS_connection_time(s)=stats(sim).num_samples_BS_connection_time(s)+1;
                   stats(sim).samples_BS_connection_time(stats(sim).num_samples_BS_connection_time(s),s)=UElist{s}{i}.end_session_time-UElist{s}{i}.start_time_serving_BS_or_relay;
                   if activate_debug_UE==1 && s==test_s && UElist{s}{i}.UEid==testUE
                        fprintf('\nTime:%d, UEid:%d STOPS. SampleTimeBS:%4.2f',time,testUE,stats(sim).samples_BS_connection_time(stats(sim).num_samples_BS_connection_time(s),s));
                   end
                end
                

                
                %Remove UE i from the list.
                UElist{s}=horzcat(UElist{s}(1:i-1),UElist{s}(i+1:num_active_UEs(s)));
                num_active_UEs(s)=num_active_UEs(s)-1;   
                %Note: the next UE to check is still the index i!!!
                %(because we have shifted the UEs in the array)
                %Then, we only increase i when we do not remove the UE.

                
            else
                i=i+1;
            end
            if i>num_active_UEs(s)
                end_process=1;
            end
        end
    end
    %CHECK SESSION STARTS
    for s=1:3
        while time_next_session_arrival(s)<=time  %By putting the while, we allow multiple arrivals in a time step.
            %Register the new UE:
            num_active_UEs(s)=num_active_UEs(s)+1;
            num_generated_UEs(s)=num_generated_UEs(s)+1;
            
            %ADD RELEVANT PARAMETERS OF THE UE, DEPENDING ON ITS TYPE
            %(Pedestrian, stationary_non_relay, stationary_relay)
            %Serving BS, Serving Relay, sp_eff, etc.
            switch s
                case 1
                    %Pedestrian
                    %Choose position inside map_valid_points_pedestrian
                    [pos_index(1),pos_index(2)]=ind2sub([sizeX,sizeY],list_valid_points_pedestrian(randi(num_valid_points_pedestrian)));
                    %Choose randomly a position inside the pixel.
                    UElist{s}{num_active_UEs(s)}.pos(1)=(pos_index(1)-1)*pixel_size+pixel_size*rand;
                    UElist{s}{num_active_UEs(s)}.pos(2)=(pos_index(2)-1)*pixel_size+pixel_size*rand;
                    
                    UElist{s}{num_active_UEs(s)}.floor=1;  %Outdoor UE, so it is in floor 1.
                    UElist{s}{num_active_UEs(s)}.end_session_time=time_next_session_arrival(s)+(-config.avg_session_duration)*log(1-rand);  
                    UElist{s}{num_active_UEs(s)}.associated_relay=0; %Only relevant for relay UEs.
                    UElist{s}{num_active_UEs(s)}.start_time=time_next_session_arrival(s);  %Start time for statistical purposes.
                    UElist{s}{num_active_UEs(s)}.speed=config.speed_pedestrian;
                    UElist{s}{num_active_UEs(s)}.indoor=0;
                    
                    if map_pedestrian_points(pos_index(1),pos_index(2))==1   %RW area
                        UElist{s}{num_active_UEs(s)}.direction=-pi+2*pi*rand;
                        UElist{s}{num_active_UEs(s)}.trajectory=0; %RW area
                    else
                        traj_number=map_numbers_traject_ped(pos_index(1),pos_index(2));
                        if traj_number>1000
                            %It is an intersection. Choose randomly between one trajectory
                            %or the other.
                            traj1=mod(traj_number,1000);
                            traj2=(traj_number-traj1)/1000;
                            if rand<0.5
                                traj_number=traj1;
                            else
                                traj_number=traj2;
                            end
                        end
                        UElist{s}{num_active_UEs(s)}.trajectory=traj_number;
                        %Choose between the direction of the trajectory or
                        %the opposite direction with prob. 0.5
                        if rand<0.5
                            UElist{s}{num_active_UEs(s)}.direction=trajectories_ped(traj_number).direction;
                        else
                            UElist{s}{num_active_UEs(s)}.direction=trajectories_ped(traj_number).direction+pi;
                        end
                        if UElist{s}{num_active_UEs(s)}.direction>pi
                            UElist{s}{num_active_UEs(s)}.direction=UElist{s}{num_active_UEs(s)}.direction-2*pi;
                        elseif UElist{s}{num_active_UEs(s)}.direction<-pi
                            UElist{s}{num_active_UEs(s)}.direction=UElist{s}{num_active_UEs(s)}.direction+2*pi;
                        end
                    end
                    
                case 2
                    %Stationary non-relay
                    
                    %Choose position inside map_valid_points_relays
                    %(since the fraction of valid points is large, we can
                    %do an iterative process)
                    valid_position=0;
                    while valid_position==0
                        %Repeat until finding a valid position
                        UElist{s}{num_active_UEs(s)}.pos(1)=scenario_x*rand;
                        UElist{s}{num_active_UEs(s)}.pos(2)=scenario_y*rand;
                        pos_index=ceil(UElist{s}{num_active_UEs(s)}.pos/pixel_size);
                        if map_valid_points_relays(pos_index(1),pos_index(2))==1
                            valid_position=1;
                        end
                    end
                    if map_buildings(pos_index(1),pos_index(2))==0
                        UElist{s}{num_active_UEs(s)}.floor=1;
                    else
                        %Indoor position: Choose randomly a floor:
                        UElist{s}{num_active_UEs(s)}.floor=randi(num_floors);
                    end
                    UElist{s}{num_active_UEs(s)}.speed=0;
                    UElist{s}{num_active_UEs(s)}.direction=0;
                    UElist{s}{num_active_UEs(s)}.trajectory=0;
                    UElist{s}{num_active_UEs(s)}.end_session_time=time_next_session_arrival(s)+(-config.avg_session_duration)*log(1-rand);
                    UElist{s}{num_active_UEs(s)}.associated_relay=0; %Only relevant for relays.
                    UElist{s}{num_active_UEs(s)}.start_time=time_next_session_arrival(s);  %Start time for statistical purposes.
                    UElist{s}{num_active_UEs(s)}.indoor=map_indoor_points(pos_index(1),pos_index(2));

            
                case 3
                    %Stationary relay
                    %Choose one of the inactive relays and associate the
                    %new UE with it.
                    list_available_relays=find(list_active_relays==0);
                    num_available_relays=size(list_available_relays,1);
                    if num_available_relays>0
                        num=list_available_relays(randi(num_available_relays));
                        list_active_relays(num)=1; 
                    
                        UElist{s}{num_active_UEs(s)}.pos=relay_list(num).pos;
                        UElist{s}{num_active_UEs(s)}.floor=relay_list(num).floor;
                        UElist{s}{num_active_UEs(s)}.speed=0;
                        UElist{s}{num_active_UEs(s)}.direction=0;
                        UElist{s}{num_active_UEs(s)}.trajectory=0;
                        UElist{s}{num_active_UEs(s)}.associated_relay=num; %This is only needed to control the activity of the relays.
                        UElist{s}{num_active_UEs(s)}.end_session_time=time_next_session_arrival(s)+(-config.avg_session_duration)*log(1-rand);
                        UElist{s}{num_active_UEs(s)}.start_time=time_next_session_arrival(s);  %Start time for statistical purposes.
                        
                        pos_index=ceil(UElist{s}{num_active_UEs(s)}.pos/pixel_size);
                        UElist{s}{num_active_UEs(s)}.indoor=map_indoor_points(pos_index(1),pos_index(2));
                    else
                        %No more relays available, so we don't activate any
                        %UE.
                        num_active_UEs(s)=num_active_UEs(s)-1;
                    end
                    
                
            end
            %Initialise serving BS and serving relay.
            UElist{s}{num_active_UEs(s)}.serving_BS=-1;
            UElist{s}{num_active_UEs(s)}.serving_relay=-1;
            UElist{s}{num_active_UEs(s)}.start_time_serving_BS_or_relay=time_next_session_arrival(s);
            
            UElist{s}{num_active_UEs(s)}.UEid=num_generated_UEs(s); %The UEid is used to track the behavior of a given UE. 
            
            if activate_debug_UE==1 && s==test_s && UElist{s}{num_active_UEs(s)}.UEid==testUE
                fprintf('\nTime:%d, UEid:%d STARTS',time,testUE);
            end

            %Schedule the next arrival
            time_next_session_arrival(s)=time_next_session_arrival(s)+(-1/config.lambda(s))*log(1-rand);
            %NOTE: We schedule not after "time" but after "time_next_session_arrival" allowing multiple sessions in a
            %time step.
    
        end
    end
    %Assess connectivity and radio conditions:
    for s=1:3
        for i=1:num_active_UEs(s)
            pos_index=ceil(UElist{s}{i}.pos/pixel_size);
            %Store the last BS and relay, in order to identify changes.
            last_BS=UElist{s}{i}.serving_BS;
            last_relay=UElist{s}{i}.serving_relay;
            if UElist{s}{i}.serving_relay>0
                connected_to_relay=1;
            else
                connected_to_relay=0;
            end
            
            if s==3
               %Stationary Relay: same conditions as relay UElist{s}{i}.associated_relay
               UElist{s}{i}.speff=relay_list(UElist{s}{i}.associated_relay).speff_BS_relay;
               UElist{s}{i}.serving_BS=relay_list(UElist{s}{i}.associated_relay).serving_BS; 
               UElist{s}{i}.serving_relay=relay_list(UElist{s}{i}.associated_relay).serving_relay;
               
            else
                %Pedestrian or stationary non-relay:
                UElist{s}{i}.speff=speff_total_relay{UElist{s}{i}.floor}(pos_index(1),pos_index(2));
                if serving_BS_relay_total{UElist{s}{i}.floor}(pos_index(1),pos_index(2))<10000
                    UElist{s}{i}.serving_BS=serving_BS_relay_total{UElist{s}{i}.floor}(pos_index(1),pos_index(2));
                    UElist{s}{i}.serving_relay=0;
                else
                    UElist{s}{i}.serving_relay=serving_BS_relay_total{UElist{s}{i}.floor}(pos_index(1),pos_index(2))-10000;
                    UElist{s}{i}.serving_BS=relay_list(UElist{s}{i}.serving_relay).serving_BS;  %BSid that serves the relay.
                    UElist{s}{i}.speff_BS=speff_total_BS{UElist{s}{i}.floor}(pos_index(1),pos_index(2));
                    UElist{s}{i}.speff_improvement=UElist{s}{i}.speff-UElist{s}{i}.speff_BS; % The spectral efficiency improvement is expressed in absolute value.                 

%                   %If the spectral efficiency improvement was computed by means of a fraction:
%                     if UElist{s}{i}.speff_BS~=0
%                         UElist{s}{i}.speff_improvement=UElist{s}{i}.speff/UElist{s}{i}.speff_BS;   
%                     else
%                         UElist{s}{i}.speff_improvement=UElist{s}{i}.speff/0.01; %This causes distortion
%                     end
                end
            end
            
            %Check if the serving BS or relay have changed and measure the
            %connection time:
            if (last_BS ~= -1) && (last_relay ~= -1)   %If both are -1 this means that it is the start of a connection.
               if connected_to_relay==1
                  %It was connected to a relay, see if now it is connected to a different relay or to a BS:
                  if (UElist{s}{i}.serving_relay ~= last_relay)  %If now it is connected to a BS, serving_relay=0, so it will be different than last_relay
                     %CHANGE =HANDOVER
                     stats(sim).num_samples_relay_connection_time(s)=stats(sim).num_samples_relay_connection_time(s)+1;
                     stats(sim).samples_relay_connection_time(stats(sim).num_samples_relay_connection_time(s),s)=time-UElist{s}{i}.start_time_serving_BS_or_relay;
                     UElist{s}{i}.start_time_serving_BS_or_relay=time;
                     if activate_debug_UE==1 && s==test_s && UElist{s}{i}.UEid==testUE
                        fprintf('\nTime:%d, UE%d, CHANGE! (WAS IN RELAY). SampleTime:%4.2f',time,testUE,stats(sim).samples_relay_connection_time(stats(sim).num_samples_relay_connection_time(s),s)); 
                     end
                  end    
               else
                   %It was connected to a BS, see if now it is connected to
                   %a relay or to a different BS. 
                    if (UElist{s}{i}.serving_relay > 0) || (UElist{s}{i}.serving_BS ~= last_BS)
                        %CHANGE=HANDOVER
                        stats(sim).num_samples_BS_connection_time(s)=stats(sim).num_samples_BS_connection_time(s)+1;
                        stats(sim).samples_BS_connection_time(stats(sim).num_samples_BS_connection_time(s),s)=time-UElist{s}{i}.start_time_serving_BS_or_relay;
                        UElist{s}{i}.start_time_serving_BS_or_relay=time;
                        if activate_debug_UE==1 && s==test_s && UElist{s}{i}.UEid==testUE
                            fprintf('\nTime:%d, UE%d, CHANGE! (WAS IN BS). SampleTime:%4.2f',time,testUE,stats(sim).samples_BS_connection_time(stats(sim).num_samples_BS_connection_time(s),s)); 
                        end
                    end
               end 
            end
            if activate_debug_UE==1 && s==test_s && UElist{s}{i}.UEid==testUE
                fprintf('\nTime:%d, UE%d, BS:%d, Relay:%d',time,testUE,UElist{s}{i}.serving_BS,UElist{s}{i}.serving_relay); 
                fprintf(' SpEff:%4.2f, PosX:%4.2f, PosY:%4.2f, Floor:%d',UElist{s}{i}.speff,UElist{s}{i}.pos(1),UElist{s}{i}.pos(2),UElist{s}{i}.floor);
            end
            
            
            %Take speff sample for UE i
            stats(sim).num_samples_speff(s)=stats(sim).num_samples_speff(s)+1;
            stats(sim).samples_speff(stats(sim).num_samples_speff(s),s)=UElist{s}{i}.speff;
            if UElist{s}{i}.serving_relay>0
                stats(sim).num_samples_connect_relay(s)=stats(sim).num_samples_connect_relay(s)+1;
            end
            if s==2 || s==3 
                if UElist{s}{i}.indoor==1
                    %Take statistic of indoor user. (note: s=2 is stored in column 1 and s=3 is stored in
                    %column 2, so the column is s-1).
                    stats(sim).num_samples_speff_indoor(s-1)=stats(sim).num_samples_speff_indoor(s-1)+1; 
                    stats(sim).samples_speff_indoor(stats(sim).num_samples_speff_indoor(s-1),s-1)=UElist{s}{i}.speff;
                    if UElist{s}{i}.serving_relay>0
                        stats(sim).num_samples_connect_relay_indoor(s-1)=stats(sim).num_samples_connect_relay_indoor(s-1)+1;
                    end
                    
                else
                    %Take statistic of outdoor user. (note: s=2 is stored in column 1 and s=3 is stored in
                    %column 2, so the column is s-1).
                    stats(sim).num_samples_speff_outdoor(s-1)=stats(sim).num_samples_speff_outdoor(s-1)+1; 
                    stats(sim).samples_speff_outdoor(stats(sim).num_samples_speff_outdoor(s-1),s-1)=UElist{s}{i}.speff;
                    if UElist{s}{i}.serving_relay>0
                        stats(sim).num_samples_connect_relay_outdoor(s-1)=stats(sim).num_samples_connect_relay_outdoor(s-1)+1;
                    end
                end
            end
        end
    
    end
    
     % Count the number of users that are connected to each relay and BS at
     % each time instant (taking into account the user category and if they
     % are indoor or outdoor in the case of stationary UEs). For BSs
     % distinguish between direct users and indirect users (connected
     % through a relay).
    
     for s = 1:3
        for i=1:num_active_UEs(s)
              
            if UElist{s}{i}.serving_relay > 0 %If the UE is being served by a relay
                this_relay=UElist{s}{i}.serving_relay;
                this_BS = UElist{s}{i}.serving_BS;
                this_BS_str = ['BS_', num2str(this_BS)];
                this_relay_str = ['relay_', num2str(this_relay)];
                if(UElist{s}{i}.indoor == 1)
                     UEs_relays{this_relay}.users{s}.indoor(1,time) = UEs_relays{this_relay}.users{s}.indoor(1,time) + 1;
                     UEs_BS.(this_BS_str){s}.indirect.indoor(1,time) = UEs_BS.(this_BS_str){s}.indirect.indoor(1,time) + 1;
                     if (UElist{s}{i}.speff_BS < config.speff_limit_outage && UElist{s}{i}.speff >= config.speff_limit_outage)
                         UEs_relays{this_relay}.users_outage{s}.indoor(1,time) = UEs_relays{this_relay}.users_outage{s}.indoor(1,time) + 1;
                     end
                else
                    UEs_relays{this_relay}.users{s}.outdoor(1,time) = UEs_relays{this_relay}.users{s}.outdoor(1,time) + 1;
                    UEs_BS.(this_BS_str){s}.indirect.outdoor(time) = UEs_BS.(this_BS_str){s}.indirect.outdoor(1,time) + 1; 
                    if (UElist{s}{i}.speff_BS < config.speff_limit_outage && UElist{s}{i}.speff >= config.speff_limit_outage)
                         UEs_relays{this_relay}.users_outage{s}.outdoor(1,time) = UEs_relays{this_relay}.users_outage{s}.outdoor(1,time) + 1;
                    end
                end
                UEs_relays{this_relay}.n_appearances = UEs_relays{this_relay}.n_appearances + 1;
%                 UEs_relays{this_relay}.speff_improvement(UEs_relays{this_relay}.n_appearances) = UElist{s}{i}.speff_improvement;
                UEs_relays{this_relay}.avg_speff_improvement = UEs_relays{this_relay}.avg_speff_improvement + UElist{s}{i}.speff_improvement;
            else
                this_BS = UElist{s}{i}.serving_BS;
                this_BS_str = ['BS_', num2str(this_BS)];
                if(UElist{s}{i}.indoor == 1)
                    UEs_BS.(this_BS_str){s}.direct.indoor(1,time) = UEs_BS.(this_BS_str){s}.direct.indoor(1,time) + 1;
                else
                    UEs_BS.(this_BS_str){s}.direct.outdoor(1,time) = UEs_BS.(this_BS_str){s}.direct.outdoor(1,time) + 1;
                end
            end        
                 
        end
     end
     
    
    %Collect statistics:
   for s=1:3
        stats(sim).num_active_UEs(time_index,s)=num_active_UEs(s);
    end
end
%Note: At the end of the simulation there may still be active UEs, but we
%don't measure the BS or relay connection times for them.
%%

% After a simulation has been completed, we create an histogram to save the
% information concerning the number of users connected to each BS or relay.
% The time dimension is lost, but otherwise it takes too much space in
% memory. The histogram indicates the number of time instants in which a
% relay has had a certain number of UEs connected (ranging from 1 to 20).

% -- Initialize histogram --
hist_size = 20; % We consider a maximum of 20 UEs connected to a relay at a given time instant
for n_relay = 1:config.num_relays
    for s=1:3
        UEs_relays{n_relay}.users{s}.indoor_hist = zeros(1,hist_size);
        UEs_relays{n_relay}.users{s}.outdoor_hist = zeros(1,hist_size);
        UEs_relays{n_relay}.users_outage{s}.indoor_hist = zeros(1,hist_size);
        UEs_relays{n_relay}.users_outage{s}.outdoor_hist = zeros(1,hist_size);
    end
end

% -- Fill histogram --
for n_relay = 1:config.num_relays
    for s=1:3
        for i = 1:hist_size
            UEs_relays{n_relay}.users{s}.indoor_hist(1,i) = length(find(UEs_relays{n_relay}.users{s}.indoor == i));
            UEs_relays{n_relay}.users{s}.outdoor_hist(1,i) = length(find(UEs_relays{n_relay}.users{s}.outdoor == i));
            UEs_relays{n_relay}.users_outage{s}.indoor_hist(1,i) = length(find(UEs_relays{n_relay}.users_outage{s}.indoor == i));
            UEs_relays{n_relay}.users_outage{s}.outdoor_hist(1,i) = length(find(UEs_relays{n_relay}.users_outage{s}.outdoor == i));
        end
    end
end

%Remove time axis data
for n_relay = 1:config.num_relays
    for s=1:3
        UEs_relays{n_relay}.users{s}= rmfield(UEs_relays{n_relay}.users{s}, 'indoor');
        UEs_relays{n_relay}.users{s}= rmfield(UEs_relays{n_relay}.users{s}, 'outdoor');
        UEs_relays{n_relay}.users_outage{s}= rmfield(UEs_relays{n_relay}.users_outage{s}, 'indoor');
        UEs_relays{n_relay}.users_outage{s}= rmfield(UEs_relays{n_relay}.users_outage{s}, 'outdoor');
    end
end

%% MEASURE STATISTICS (PER SIMULATION)
stats(sim).avg_active_UEs=mean(stats(sim).num_active_UEs);

size_samples_speff=max([stats(sim).num_samples_speff(1),stats(sim).num_samples_speff(2),stats(sim).num_samples_speff(3)]);
size_samples_BS_connection_time=max([stats(sim).num_samples_BS_connection_time(1),stats(sim).num_samples_BS_connection_time(2),stats(sim).num_samples_BS_connection_time(3)]);
size_samples_relay_connection_time=max([stats(sim).num_samples_relay_connection_time(1),stats(sim).num_samples_relay_connection_time(2),stats(sim).num_samples_relay_connection_time(3)]);
stats(sim).samples_speff=stats(sim).samples_speff(1:size_samples_speff,:);
stats(sim).samples_BS_connection_time=stats(sim).samples_BS_connection_time(1:size_samples_BS_connection_time,:);
stats(sim).samples_relay_connection_time=stats(sim).samples_relay_connection_time(1:size_samples_relay_connection_time,:);
for s=1:3
    stats(sim).avg_session_duration(s)=stats(sim).avg_session_duration(s)/stats(sim).num_sessions_released(s);
    
    %SpEff statistics: we don't have to consider the samples equal to 0 due
    %to the initial preallocation.
    stats(sim).speff_avg(s)=mean(stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s));
    stats(sim).speff_perc5(s)=prctile(stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s),5);
    stats(sim).speff_perc95(s)=prctile(stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s),95);
    
    auxstats=stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s);
    stats(sim).outage_prob(s)=size(auxstats(auxstats<=config.speff_limit_outage),1)/size(auxstats,1);
    stats(sim).cdf_speff{s}=generate_CDF(auxstats);
       
    stats(sim).BS_connection_time_avg(s)=mean(stats(sim).samples_BS_connection_time(1:stats(sim).num_samples_BS_connection_time(s),s));
    stats(sim).BS_connection_time_perc5(s)=prctile(stats(sim).samples_BS_connection_time(1:stats(sim).num_samples_BS_connection_time(s),s),5);
    stats(sim).BS_connection_time_perc95(s)=prctile(stats(sim).samples_BS_connection_time(1:stats(sim).num_samples_BS_connection_time(s),s),95);
    stats(sim).cdf_BS_connection_time{s}=generate_CDF(stats(sim).samples_BS_connection_time(1:stats(sim).num_samples_BS_connection_time(s),s));
    
    stats(sim).relay_connection_time_avg(s)=mean(stats(sim).samples_relay_connection_time(1:stats(sim).num_samples_relay_connection_time(s),s));
    stats(sim).relay_connection_time_perc5(s)=prctile(stats(sim).samples_relay_connection_time(1:stats(sim).num_samples_relay_connection_time(s),s),5);
    stats(sim).relay_connection_time_perc95(s)=prctile(stats(sim).samples_relay_connection_time(1:stats(sim).num_samples_relay_connection_time(s),s),95);
    stats(sim).cdf_relay_connection_time{s}=generate_CDF(stats(sim).samples_relay_connection_time(1:stats(sim).num_samples_relay_connection_time(s),s));
    
    stats(sim).prob_connect_relay(s)=stats(sim).num_samples_connect_relay(s)/stats(sim).num_samples_speff(s);
end

for s=1:2  %Note that for indoor/outdoor stats we only consider two categories of UEs (the stationary UEs non-relay and relay)
    stats(sim).speff_indoor_avg(s)=mean(stats(sim).samples_speff_indoor(1:stats(sim).num_samples_speff_indoor(s),s));
    stats(sim).speff_indoor_perc5(s)=prctile(stats(sim).samples_speff_indoor(1:stats(sim).num_samples_speff_indoor(s),s),5);
    stats(sim).speff_indoor_perc95(s)=prctile(stats(sim).samples_speff_indoor(1:stats(sim).num_samples_speff_indoor(s),s),95);
    
    auxstats=stats(sim).samples_speff_indoor(1:stats(sim).num_samples_speff_indoor(s),s);
    stats(sim).outage_prob_indoor(s)=size(auxstats(auxstats<=config.speff_limit_outage),1)/size(auxstats,1);
    stats(sim).cdf_speff_indoor{s}=generate_CDF(auxstats);
    
    stats(sim).speff_outdoor_avg(s)=mean(stats(sim).samples_speff_outdoor(1:stats(sim).num_samples_speff_outdoor(s),s));
    stats(sim).speff_outdoor_perc5(s)=prctile(stats(sim).samples_speff_outdoor(1:stats(sim).num_samples_speff_outdoor(s),s),5);
    stats(sim).speff_outdoor_perc95(s)=prctile(stats(sim).samples_speff_outdoor(1:stats(sim).num_samples_speff_outdoor(s),s),95);
    
    auxstats=stats(sim).samples_speff_outdoor(1:stats(sim).num_samples_speff_outdoor(s),s);
    stats(sim).outage_prob_outdoor(s)=size(auxstats(auxstats<=config.speff_limit_outage),1)/size(auxstats,1);
    stats(sim).cdf_speff_outdoor{s}=generate_CDF(auxstats);
    
    stats(sim).prob_connect_relay_indoor(s)=stats(sim).num_samples_connect_relay_indoor(s)/stats(sim).num_samples_speff_indoor(s);
    stats(sim).prob_connect_relay_outdoor(s)=stats(sim).num_samples_connect_relay_outdoor(s)/stats(sim).num_samples_speff_outdoor(s);
end

%Statistics aggregated for all categories:
vector_aux=stats(sim).samples_speff(1:stats(sim).num_samples_speff(1),1);
for s=2:3
    vector_aux=[vector_aux;stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s)];
end
stats(sim).speff_avg_total=mean(vector_aux);
stats(sim).speff_perc5_total=prctile(vector_aux,5);
stats(sim).speff_perc95_total=prctile(vector_aux,95);
stats(sim).outage_prob_total=size(vector_aux(vector_aux<=config.speff_limit_outage),1)/size(vector_aux,1);
stats(sim).cdf_speff_total=generate_CDF(vector_aux);    


stats(sim).UEs_BS = UEs_BS;
if(parameter_to_change ~= 0)
    for n_relay = 1:config.num_relays
        if UEs_relays{n_relay}.n_appearances ~= 0
            UEs_relays{n_relay}.avg_speff_improvement = (UEs_relays{n_relay}.avg_speff_improvement/UEs_relays{n_relay}.n_appearances);
        else
            UEs_relays{n_relay}.avg_speff_improvement = 0;
        end
    end
    stats(sim).UEs_relays = UEs_relays;
end
% 
end %End of simulations

%% MEASURE GLOBAL STATISTICS BY AVERAGING ALL THE SIMULATIONS
for s=1:3
    vector_aux=stats(1).num_active_UEs(:,s);
    for sim=2:config.num_simulations
        vector_aux=[vector_aux;stats(sim).num_active_UEs(:,s)];
    end
    statsglobal.avg_active_UEs(s)=mean(vector_aux);
    
    sum_num=0;
    sum_denom=0;
    for sim=1:config.num_simulations
       sum_num=sum_num+stats(sim).avg_session_duration(s)*stats(sim).num_sessions_released(s);
       sum_denom=sum_denom+stats(sim).num_sessions_released(s);
    end
    statsglobal.avg_session_duration(s)=sum_num/sum_denom;
    
    vector_aux=stats(1).samples_speff(1:stats(1).num_samples_speff(s),s);
    for sim=2:config.num_simulations
        vector_aux=[vector_aux;stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s)];
    end
    statsglobal.speff_avg(s)=mean(vector_aux);
    statsglobal.speff_perc5(s)=prctile(vector_aux,5);
    statsglobal.speff_perc95(s)=prctile(vector_aux,95);
    statsglobal.outage_prob(s)=size(vector_aux(vector_aux<=config.speff_limit_outage),1)/size(vector_aux,1);
    statsglobal.cdf_speff{s}=generate_CDF(vector_aux);
    
    vector_aux=stats(1).samples_BS_connection_time(1:stats(1).num_samples_BS_connection_time(s),s);
    for sim=2:config.num_simulations
        vector_aux=[vector_aux;stats(sim).samples_BS_connection_time(1:stats(sim).num_samples_BS_connection_time(s),s)];
    end
    statsglobal.BS_connection_time_avg(s)=mean(vector_aux);
    statsglobal.BS_connection_time_perc5(s)=prctile(vector_aux,5);
    statsglobal.BS_connection_time_perc95(s)=prctile(vector_aux,95);
    statsglobal.cdf_BS_connection_time{s}=generate_CDF(vector_aux);
    
    vector_aux=stats(1).samples_relay_connection_time(1:stats(1).num_samples_relay_connection_time(s),s);
    for sim=2:config.num_simulations
        vector_aux=[vector_aux;stats(sim).samples_relay_connection_time(1:stats(sim).num_samples_relay_connection_time(s),s)];
    end
    statsglobal.relay_connection_time_avg(s)=mean(vector_aux);
    statsglobal.relay_connection_time_perc5(s)=prctile(vector_aux,5);
    statsglobal.relay_connection_time_perc95(s)=prctile(vector_aux,95);
    statsglobal.cdf_relay_connection_time{s}=generate_CDF(vector_aux);
    
    sum_num=0;
    sum_denom=0;
    for sim=1:config.num_simulations
       sum_num=sum_num+stats(sim).num_samples_connect_relay(s);
       sum_denom=sum_denom+stats(sim).num_samples_speff(s);
    end
    statsglobal.prob_connect_relay(s)=sum_num/sum_denom;
    
    
end

for s=1:2
   %Statistics of outdoor/indoor 
   vector_aux=stats(1).samples_speff_indoor(1:stats(1).num_samples_speff_indoor(s),s);
   for sim=2:config.num_simulations
       vector_aux=[vector_aux;stats(sim).samples_speff_indoor(1:stats(sim).num_samples_speff_indoor(s),s)];
   end
   statsglobal.speff_indoor_avg(s)=mean(vector_aux);
   statsglobal.speff_indoor_perc5(s)=prctile(vector_aux,5);
   statsglobal.speff_indoor_perc95(s)=prctile(vector_aux,95);
   statsglobal.outage_prob_indoor(s)=size(vector_aux(vector_aux<=config.speff_limit_outage),1)/size(vector_aux,1);
   statsglobal.cdf_speff_indoor{s}=generate_CDF(vector_aux);
   
   vector_aux=stats(1).samples_speff_outdoor(1:stats(1).num_samples_speff_outdoor(s),s);
   for sim=2:config.num_simulations
       vector_aux=[vector_aux;stats(sim).samples_speff_outdoor(1:stats(sim).num_samples_speff_outdoor(s),s)];
   end
   statsglobal.speff_outdoor_avg(s)=mean(vector_aux);
   statsglobal.speff_outdoor_perc5(s)=prctile(vector_aux,5);
   statsglobal.speff_outdoor_perc95(s)=prctile(vector_aux,95);
   statsglobal.outage_prob_outdoor(s)=size(vector_aux(vector_aux<=config.speff_limit_outage),1)/size(vector_aux,1);
   statsglobal.cdf_speff_outdoor{s}=generate_CDF(vector_aux);
   
   sum_num=0;
   sum_denom=0;
   for sim=1:config.num_simulations
      sum_num=sum_num+stats(sim).num_samples_connect_relay_indoor(s);
      sum_denom=sum_denom+stats(sim).num_samples_speff_indoor(s);
   end
   statsglobal.prob_connect_relay_indoor(s)=sum_num/sum_denom;
   
   sum_num=0;
   sum_denom=0;
   for sim=1:config.num_simulations
      sum_num=sum_num+stats(sim).num_samples_connect_relay_outdoor(s);
      sum_denom=sum_denom+stats(sim).num_samples_speff_outdoor(s);
   end
   statsglobal.prob_connect_relay_outdoor(s)=sum_num/sum_denom;
   
end
%Statistics aggregated for all categories:
vector_aux=[];
for sim=1:config.num_simulations
   for s=1:3
       vector_aux=[vector_aux;stats(sim).samples_speff(1:stats(sim).num_samples_speff(s),s)];
   end
end   
statsglobal.speff_avg_total=mean(vector_aux);
statsglobal.speff_perc5_total=prctile(vector_aux,5);
statsglobal.speff_perc95_total=prctile(vector_aux,95);
statsglobal.outage_prob_total=size(vector_aux(vector_aux<=config.speff_limit_outage),1)/size(vector_aux,1);
statsglobal.cdf_speff_total=generate_CDF(vector_aux);    

% In order to save the file, let's remove all the vectors of samples in
% stats:
if config.store_samples_in_file==0
for sim=1:config.num_simulations
   %Don't include the samples of num_active_UEs, speff, BS_connection_time, relay_connection_time.
   %Don't include the individual CDFs per simulation, since we have the
   %global one.
   stats_to_store(sim).avg_session_duration=stats(sim).avg_session_duration;
   stats_to_store(sim).num_sessions_released=stats(sim).num_sessions_released;
   stats_to_store(sim).num_samples_speff=stats(sim).num_samples_speff;
   stats_to_store(sim).num_samples_relay_connection_time=stats(sim).num_samples_relay_connection_time;
   stats_to_store(sim).num_samples_BS_connection_time=stats(sim).num_samples_BS_connection_time;
   stats_to_store(sim).num_samples_connect_relay=stats(sim).num_samples_connect_relay;
   stats_to_store(sim).avg_active_UEs=stats(sim).avg_active_UEs;
   stats_to_store(sim).speff_avg=stats(sim).speff_avg;
   stats_to_store(sim).speff_perc5=stats(sim).speff_perc5;
   stats_to_store(sim).speff_perc95=stats(sim).speff_perc95;
   stats_to_store(sim).outage_prob=stats(sim).outage_prob;
   stats_to_store(sim).BS_connection_time_avg=stats(sim).BS_connection_time_avg;
   stats_to_store(sim).BS_connection_time_perc5=stats(sim).BS_connection_time_perc5;
   stats_to_store(sim).BS_connection_time_perc95=stats(sim).BS_connection_time_perc95;
   stats_to_store(sim).relay_connection_time_avg=stats(sim).relay_connection_time_avg;
   stats_to_store(sim).relay_connection_time_perc5=stats(sim).relay_connection_time_perc5;
   stats_to_store(sim).relay_connection_time_perc95=stats(sim).relay_connection_time_perc95;
   stats_to_store(sim).prob_connect_relay=stats(sim).prob_connect_relay;
   if config.num_relays>0
        stats_to_store(sim).relay_list=stats(sim).relay_list;
   end
   stats_to_store(sim).speff_indoor_avg=stats(sim).speff_indoor_avg;
   stats_to_store(sim).speff_indoor_perc5=stats(sim).speff_indoor_perc5;
   stats_to_store(sim).speff_indoor_perc95=stats(sim).speff_indoor_perc95;
   stats_to_store(sim).outage_prob_indoor=stats(sim).outage_prob_indoor;
   
   stats_to_store(sim).speff_outdoor_avg=stats(sim).speff_outdoor_avg;
   stats_to_store(sim).speff_outdoor_perc5=stats(sim).speff_outdoor_perc5;
   stats_to_store(sim).speff_outdoor_perc95=stats(sim).speff_outdoor_perc95;
   stats_to_store(sim).outage_prob_outdoor=stats(sim).outage_prob_outdoor; 
   
   stats_to_store(sim).prob_connect_relay_indoor=stats(sim).prob_connect_relay_indoor;
   stats_to_store(sim).prob_connect_relay_outdoor=stats(sim).prob_connect_relay_outdoor;
   
   stats_to_store(sim).speff_avg_total=stats(sim).speff_avg_total;
   stats_to_store(sim).speff_perc5_total=stats(sim).speff_perc5_total;
   stats_to_store(sim).speff_perc95_total=stats(sim).speff_perc95_total;
   stats_to_store(sim).outage_prob_total=stats(sim).outage_prob_total;
   
   if(parameter_to_change ~= 0)
    stats_to_store(sim).UEs_relays=stats(sim).UEs_relays;
   end
   stats_to_store(sim).UEs_BS=stats(sim).UEs_BS;
end

stats=stats_to_store;
clear stats_to_store;
end
% 

%save('Results.mat','-v7.3');
save(config.output_file_name,'stats','statsglobal','config', '-v7.3');
% save(config.output_file_name,'config', '-v7.3');
if exist(config.map_BS_file_name)==0
    %Store the file with the maps of the BSs if it has not been stored
    %previously.
    save(config.map_BS_file_name,'BS','SNR_total_BS','speff_total_BS','serving_BS','config');
end

end


