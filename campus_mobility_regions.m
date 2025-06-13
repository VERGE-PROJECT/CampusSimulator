%%%THIS SCRIPT COMPUTES THE MAPS OF THE DIFFERENT MOBILITY REGIONS
%%% It is called by the 'create_scenario.mat' script:

%%%Specifically, the relevant maps are the following ones:
% - map_trajectories_ped: Binary map with the trajectories of pedestrians in the
% strets
% - map_numbers_traject_ped: Map that indicates in each pixel the number of
% trajectory. If it is an intersection between traj n,m, it is encoded as
% 1000*n+m  where n>m
% - map_pedestrian_points: Binary map with the positions where pedestrians
% move with random walk.
% - map_valid_points_pedestrian: Binary map with the positions where
% pedestrians can move (trajectories + random walk areas)
% - list_valid_points_pedestrian: List with the linear indices of the valid
% points for pedestrians.
% - num_valid_points_pedestrian: Number of valid points for pedestrians.
%
% - map_valid_points_relays: Binary map with the positions where we can
% have a fixed relay 
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%PEDESTRIAN TRAJECTORIES: mobility according to direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- HORIZONTAL --
trajectories_ped(1).points=[0,3; 336,3]; %Trajectory goes from point in first row to point in second row.
trajectories_ped(1).width=6;  %It is +/- width/2 around the points
trajectories_ped(2).points=[0,30; 336,30]; 
trajectories_ped(2).width=10;  
trajectories_ped(3).points=[0,61; 336,61]; 
trajectories_ped(3).width=12;
trajectories_ped(4).points=[0,93; 336,93]; 
trajectories_ped(4).width=14;
trajectories_ped(5).points=[0,122; 336,122]; 
trajectories_ped(5).width=6;
% -- VERTICAL --
trajectories_ped(6).points=[2.5,0; 2.5,126]; 
trajectories_ped(6).width=5;
trajectories_ped(7).points=[109.5,0; 109.5,126]; 
trajectories_ped(7).width=5;
trajectories_ped(8).points=[153.5,0; 153.5,126]; 
trajectories_ped(8).width=5;
trajectories_ped(9).points=[197.5,0; 197.5,126]; 
trajectories_ped(9).width=5;
trajectories_ped(10).points=[242.5,0; 242.5,126]; 
trajectories_ped(10).width=5;
trajectories_ped(11).points=[285.5,0; 285.5,126]; 
trajectories_ped(11).width=5;
trajectories_ped(12).points=[331,0; 331,126];
trajectories_ped(12).width=8;

num_trajectories_ped=size(trajectories_ped,2);

for n=1:num_trajectories_ped
    trajectories_ped(n).direction=atan2(trajectories_ped(n).points(2,2)-trajectories_ped(n).points(1,2),trajectories_ped(n).points(2,1)-trajectories_ped(n).points(1,1));
end

[map_trajectories_ped,map_numbers_traject_ped] = fill_trajectory_map(sizeX,sizeY,pixel_size,map_indoor_points,trajectories_ped);

map_buildings_with_traj=map_indoor_points+2*map_trajectories_ped;
%figure; imshow((map_trajectories_ped'),[0,1]);colormap(jet);colorbar;title('Trajectories Pedestrians');
%figure; imshow((map_numbers_traject_ped'),[0,num_trajectories_ped]);colormap(jet);colorbar;title('Trajectories Pedestrians');
%figure; imshow((map_buildings_with_traj'),[0,2]);colormap(jet);colorbar;title('Buildings and Trajectories');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%PEDESTRIAN AREAS: mobility according to random walk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pedestrian_areas(1).vertices=[90,67; 129,67; 129,86; 90,86];
pedestrian_areas(2).vertices=[178,67; 217,67; 217,86; 178,86];
pedestrian_areas(3).vertices=[266,67; 306,67; 306,86; 266,86];
pedestrian_areas(4).vertices=[5,35; 5,55; 18,55; 18,35];
pedestrian_areas(5).vertices=[5,67; 50,67; 50,86; 5,86];
pedestrian_areas(6).vertices=[5,100; 68,100; 68,119; 5,119];
% 
[map_pedestrian_points,~]=find_indoor_points(sizeX,sizeY,pixel_size,pedestrian_areas);
map_pedestrian_points=map_pedestrian_points.*(1-map_indoor_points); %This to force that all pedestrian points are outdoor.
%figure; imshow((map_pedestrian_points'),[0,1]);colormap(jet);colorbar;title('Pedestrian areas');

map_buildings_with_traj_and_ped_areas=map_buildings_with_traj+3*map_pedestrian_points;
%figure; imshow((map_buildings_with_traj_and_ped_areas'),[0,3]);colormap(jet);colorbar;title('Buildings,Traject. and pedest. areas');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%MAPS OF VALID POINTS FOR MOBILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Map of positions where the pedestrians can move: the trajectories and the
%random walk areas (no need to exclude the indoor points because they are
%not included in these maps).
map_valid_points_pedestrian=min(1,map_trajectories_ped+map_pedestrian_points);
list_valid_points_pedestrian=find(map_valid_points_pedestrian==1);
num_valid_points_pedestrian=size(list_valid_points_pedestrian,1);

%Map of positions where the relays can be: all the points (indoor/outdoor)
map_valid_points_relays=ones(sizeX,sizeY);

%figure; imshow((map_valid_points_pedestrian'),[0,1]);colormap(jet);colorbar;title('Map Valid Points Pedestrians');
%figure; imshow((map_valid_points_relays'),[0,1]);colormap(jet);colorbar;title('Map Valid Points Relays');


