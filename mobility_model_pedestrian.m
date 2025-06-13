function [new_pos,new_direction,new_traj] = mobility_model_pedestrian(pos,direction,traj,speed,time_step,prob_change_RW,angle_change_RW,prob_change_inters,sizeX,sizeY,pixel_size,map_valid_points,map_RW,trajectories,map_trajectories)
%MOBILITY_MODEL_PEDESTRIAN:  Executes the movement of one pedestrian UE in
%one time_step. Input parameters:
%   -pos: Initial position: vector with the x,y component in m 
%   -direction: in rad (between -pi and pi)
%   -traj: Identifier of the current trajectory (0 if it is in a RW area)
%   -speed: in m/s
%   -time_step: duration of the time step in s.
%   -prob_change_RW: probability of changing direction in the time step
%   with random walk
%   - angle_change_RW: angle of direction change (between +angle_change_RW
%   and -angle_change_RW) measured in rad. Used by Random Walk.
%   - prob_change_inters: probability of changing direction when reaching
%   an intersection in a street.
%   - sizeX,sizeY,pixel_size: parameters of the maps
%   - map_valid_points: Map with the valid points for a pedestrian (marked 1)
%   - map_RW: Map where the UE moves using random walk (pedestrian areas)
%   - trajectories: List of trajectories characterized by points, widths and
%   directions
%   - map_trajectories: Map with the trajectory associated to each pixel.
%
%   OUTPUTS:
%   -new_pos: New position
%   -new_direction: New direction
%   -new_traj: Identifier of the new trajectory.

pos_index=ceil(pos/pixel_size);
if map_valid_points(pos_index(1),pos_index(2))==0
    fprintf('\nERROR: The UE is in an invalid position!!');
    %Do nothing
    new_pos=pos;
    new_direction=direction;
    new_traj=traj;
    
elseif map_RW(pos_index(1),pos_index(2))==1
    %Random walk area
    new_pos=pos+speed*time_step*[cos(direction),sin(direction)];
    new_pos_index=ceil(new_pos/pixel_size);
    
    %Posibilities to check:
    %1) The UE goes out of the scenario
    %2) The UE enters an invalid point
    %In any of the cases it goes back in the opposite direction: new_dir=dir-pi
    if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
        %Note: the last expression of the if is not evaluated if any of the
        %previous ones is true, so it is not possible to evaluate the map
        %in a point outside the indices.
        
        new_direction=direction-pi;
        
        %Update the position:
        new_pos=pos+speed*time_step*[cos(new_direction),sin(new_direction)];
        new_pos_index=ceil(new_pos/pixel_size);
        
        %In some cases the new position can still be invalid. If so, make a
        %movement in either new_direction+pi/2 or new_direction-pi/2.
        if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
            new_pos=pos+speed*time_step*[cos(new_direction+0.5*pi),sin(new_direction+0.5*pi)];
            new_pos_index=ceil(new_pos/pixel_size);
            if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                new_pos=pos+speed*time_step*[cos(new_direction-0.5*pi),sin(new_direction-0.5*pi)];
                new_pos_index=ceil(new_pos/pixel_size);
                %In the very specific case that the position is still invalid, simply try with the four
                %adjacent pixels to the original position, and if not, keep the UE at the original position.
                if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                    new_pos_index=pos_index+[1,0];
                    new_pos=pos+pixel_size*[1,0];
                    if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                        new_pos_index=pos_index+[-1,0];
                        new_pos=pos+pixel_size*[-1,0];
                        if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                            new_pos_index=pos_index+[0,1];
                            new_pos=pos+pixel_size*[0,1];
                            if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                                new_pos_index=pos_index+[0,-1];
                                new_pos=pos+pixel_size*[0,-1];
                                if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                                    new_pos_index=pos_index;
                                    new_pos=pos;
                                end
                            end
                        end
                    end
                end
                
                
            end
        end
        
        
    else
       new_direction=direction; %For the time being, let's keep the old direction, perhaps it will be changed later.
    end
    
    %Direction update:
    if map_RW(new_pos_index(1),new_pos_index(2))==1
        %The UE remains in the random walk area.
        %Check if there is a change in direction
        if rand<prob_change_RW
            %Change in direction in a margin
            %new_direction+unif(-angle_change_RW, angle_change_RW)
            new_direction=new_direction-angle_change_RW + 2*angle_change_RW*rand;
        end
        new_traj=0;
    else
        %It goes out of the RW area, so necessarily it must be in a
        %pedestrian trajectory (otherwise there will be an error...)
        traj_number=map_trajectories(new_pos_index(1),new_pos_index(2));
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
        new_traj=traj_number;
        %With prob 0.5 we choose the direction of the trajectory and
        %with prob 0.5 the opposite direction
        if rand<0.5
            new_direction=trajectories(traj_number).direction;
        else
        	new_direction=trajectories(traj_number).direction+pi;
        end
    end
else
    %Movement in a street.
    new_traj=traj;  %Initialise trajectory to the original value.
    %Update in the direction
    new_pos=pos+speed*time_step*[cos(direction),sin(direction)];
    new_pos_index=ceil(new_pos/pixel_size);
    
     
    %Check if we exceed the limits of the scenario or enter an invalid
    %zone. If so, just go back.
    if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
        new_direction=direction-pi;
        
        %Update the position:
        new_pos=pos+speed*time_step*[cos(new_direction),sin(new_direction)];
        new_pos_index=ceil(new_pos/pixel_size);
        
        %In some cases the new position can still be invalid. If so, make a
        %movement in either new_direction+pi/2 or new_direction-pi/2.
        %However, we maintain the direction of the street (i.e.
        %direction-pi)
        if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
            new_pos=pos+speed*time_step*[cos(new_direction+0.5*pi),sin(new_direction+0.5*pi)];
            new_pos_index=ceil(new_pos/pixel_size);
            if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                new_pos=pos+speed*time_step*[cos(new_direction-0.5*pi),sin(new_direction-0.5*pi)];
                new_pos_index=ceil(new_pos/pixel_size);
                
                %In the very specific case that the position is still invalid, simply try with the four
                %adjacent pixels to the original position, and if not, keep the UE at the original position.
                if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                    new_pos_index=pos_index+[1,0];
                    new_pos=pos+pixel_size*[1,0];
                    if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                        new_pos_index=pos_index+[-1,0];
                        new_pos=pos+pixel_size*[-1,0];
                        if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                            new_pos_index=pos_index+[0,1];
                            new_pos=pos+pixel_size*[0,1];
                            if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                                new_pos_index=pos_index+[0,-1];
                                new_pos=pos+pixel_size*[0,-1];
                                if new_pos_index(1)<1 || new_pos_index(1)>sizeX || new_pos_index(2)<1 || new_pos_index(2)>sizeY || map_valid_points(new_pos_index(1),new_pos_index(2))==0
                                    new_pos_index=pos_index;
                                    new_pos=pos;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        new_direction=direction; %For the time being, let's keep the old direction
    end
    
    %Check if it is an intersection to change direction (note that if the
    %new point is in a RW area the new trajectory will not really have
    %effect because in the next step if will move according to RW model)
    %if map_RW(pos_index(1),pos_index(2))==0
        traj_number=map_trajectories(new_pos_index(1),new_pos_index(2));
        if traj_number>1000
           %It is an intersection. 
           if rand<prob_change_inters
              %Change trajectory and choose randomly one of the two possible directions
              traj1=mod(traj_number,1000);
              traj2=(traj_number-traj1)/1000;
              %if abs(direction-trajectories(traj1).direction)<1E-6 || abs(abs(direction-trajectories(traj1).direction)/pi-1)<1E-6
              %if direction==trajectories(traj1).direction || direction==trajectories(traj1).direction+pi || direction==trajectories(traj1).direction-pi
                  %Note that this check is performed with the initial
                  %direction, because the new direction may have changed if
                  %turning back due to exceeding limits, and if so it may
                  %happen that the new one is not between pi and -pi
              if traj==traj1    
                  new_traj=traj2;  %It is in traj1 so now changes to traj2
              else
                  new_traj=traj1; %It is in traj2 so now changes to traj1
              end
              if rand<0.5
                 new_direction=trajectories(new_traj).direction;
              else
                 new_direction=trajectories(new_traj).direction+pi; 
              end
           end
        end
    %else
        %New position is in RW area, so trajectory is now 0.
    %    new_traj=0;
    %end
    
    
    
end

%Make sure the new direction is between -pi and pi
if new_direction>pi
    new_direction=new_direction-2*pi;
elseif new_direction<-pi
    new_direction=new_direction+2*pi;
end
 
end

