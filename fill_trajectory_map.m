function [map_trajectories,map_numbers] = fill_trajectory_map(sizeX,sizeY,pixel_size,map_indoor_points,trajectories)
%FILL_TRAJECTORY_MAP Summary of this function goes here
%   This function creates two maps:
%   - map_trajectories: Value 1 if the pixel corresponds to a trajectory
%   - map_numbers: Number of the trajectory of the pixel if there is one
%   trajectory. If it is an intersection of 2 traj. of numbers n1,n2 they are encoded as
%   n1*1000+n2
%   
%   Inputs:
%   -sizeX,sizeY,pixel_size: parameters of the scenario.
%   -map_indoor_points: Map with the pixels that are indoor. To avoid that
%   they are selected as trajectory points.
%   -trajectories: struct with the trajectories. It includes the field
%   ".points" with the edge points of the trajectory

num_trajectories=size(trajectories,2);
map_trajectories=zeros(sizeX,sizeY);
map_numbers=zeros(sizeX,sizeY);

for n=1:num_trajectories
    step=pixel_size;
    trajectory_length=norm(trajectories(n).points(2,:)-trajectories(n).points(1,:));
    covered_length=0;
    pos=trajectories(n).points(1,:); %Start point
    while covered_length<trajectory_length
        %Fill the pixels in the width direction (direction-pi/2)
        %corresponding to a line orthogonal to the direction and centered
        %at pos
        pos_width=pos+0.5*trajectories(n).width*[cos(trajectories(n).direction+pi/2),sin(trajectories(n).direction+pi/2)];
        covered_width=0;
        while covered_width<trajectories(n).width
            %pos_index=[max(1,min(ceil(pos_width(1)/pixel_size),sizeX)),max(1,min(ceil(pos_width(2)/pixel_size),sizeY))];%This verifies it is in the limits of the scenario
            pos_index=ceil(pos_width/pixel_size);
            if pos_index(1)>=1 && pos_index(1)<=sizeX && pos_index(2)>=1 && pos_index(2)<=sizeY && map_indoor_points(pos_index(1),pos_index(2))==0
                %To make sure that only points inside the scenario and not
                %indoor are selected.
            %if map_indoor_points(pos_index(1),pos_index(2))==0   %To make sure that no indoor point is selected as part of the trajectory.
                map_trajectories(pos_index(1),pos_index(2))=1; 
                if map_numbers(pos_index(1),pos_index(2))==0
                    map_numbers(pos_index(1),pos_index(2))=n;
                elseif map_numbers(pos_index(1),pos_index(2))<n
                    %This point is an intersection.
                    %Encode the intersection between n1 and n2 (n2>n1) as:
                    %n2*1000+n1
                    map_numbers(pos_index(1),pos_index(2))=n*1000+map_numbers(pos_index(1),pos_index(2));
                    
                    %Note: if the stored number is >=n it means that we
                    %have previously marked this point as either belonging
                    %to trajectory n (only this itrajectory if equal to n,
                    %or as an intersection if higher than n).
                end
            end
            pos_width=pos_width+step*0.5*[cos(trajectories(n).direction-pi/2),sin(trajectories(n).direction-pi/2)];
            covered_width=covered_width+step*0.5;
        end
        pos=pos+step*[cos(trajectories(n).direction),sin(trajectories(n).direction)];
        covered_length=covered_length+step;
    end 
end
end

