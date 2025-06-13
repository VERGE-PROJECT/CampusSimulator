function map_indoor_distances= find_indoor_distances_relay(sizeX,sizeY,pixel_size,posBS,buildings,map_buildings,dmax)
%FIND_INDOOR_DISTANCES Summary of this function goes here
%   This function computes the 2D indoor distance between a relay and another point.
%   If the relay is indoor the distance is from the relay to the wall of
%   its own building.
%   If the relay is outdoor the distance is from an indoor point to the
%   wall of the building of this point
%   Computations are only done for points that are at distance lower than
%   dmax.
%Inputs:
%- Scenario parameters (sizeX,sizeY) in pixels, pixel_size in m
%- posBS: position of the BS
%- buildings: struct with the vertices of each building
%- map_buildings: Map that indicates the building where each pixel belongs
%(0 if it is an outdoor point)
%- dmax: maximum distance for the computation



pixels_max=ceil(dmax/pixel_size);
posBS_index=ceil(posBS/pixel_size);

building_relay=map_buildings(posBS_index(1),posBS_index(2));


if building_relay==0
    %Outdoor relay:
    %Map will be 0 for outdoor points and distance for buildings.
    %Initial map: distance=infinite (1E6) inside buildings.
    map_indoor_distances=min(1,map_buildings)*1E6; 
    

    %Check only a square centered in the position of the relay and distance
    %2*pixels_max. All the other indoor points will maintain infinite distance.
    for i=max(1,posBS_index(1)-pixels_max):min(sizeX,posBS_index(1)+pixels_max)
        for j=max(1,posBS_index(2)-pixels_max):min(sizeY,posBS_index(2)+pixels_max)
            pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
            
            b=map_buildings(i,j);
            if b>0 && norm(pos-posBS)<=dmax
                %Indoor point
                %Method:
                %- Check intersections of line pos-posBS with all the segments around the
                %building.
                %- When there is an intersection, compute the distance between
                %pos and the intersection point.
                %- For buildings with "normal" shapes there will be 1
                %intersection. In some cases of "strange" shapes there could be
                %more than one. In this case, we take the minimum distance.
           
                num_vertices=size(buildings(b).vertices,1);
                for v=1:num_vertices
                    next_v=mod(v,num_vertices)+1;
                    [intersection,point]=intersect_lines(posBS(1),posBS(2),pos(1),pos(2),buildings(b).vertices(v,1),buildings(b).vertices(v,2),buildings(b).vertices(next_v,1),buildings(b).vertices(next_v,2));
                    if intersection==1
                        %There is an intersection
                        map_indoor_distances(i,j)=min(norm(pos-point),map_indoor_distances(i,j));
                    end
                end
             end
        end
    end
else
    %Indoor relay:
    %Map will be the distance from the relay to the intersection between the wall and the line between a point and the relay
    %Only computed for outdoor points. For indoor points it is 0.
    %Initial map distance=infinite (1E6) in outdoor points.
    map_indoor_distances=(1-min(1,map_buildings))*1E6; 
    %Check only a square centered in the position of the relay and distance
    %2*pixels_max. All the other indoor points will maintain infinite distance.
    for i=max(1,posBS_index(1)-pixels_max):min(sizeX,posBS_index(1)+pixels_max)
        for j=max(1,posBS_index(2)-pixels_max):min(sizeY,posBS_index(2)+pixels_max)
           pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
           if map_buildings(i,j)==0 && norm(pos-posBS)<=dmax
               %Outdoor point
               %Check intersections with the walls of building_relay
               num_vertices=size(buildings(building_relay).vertices,1);
               for v=1:num_vertices
                    next_v=mod(v,num_vertices)+1;
                    [intersection,point]=intersect_lines(posBS(1),posBS(2),pos(1),pos(2),buildings(building_relay).vertices(v,1),buildings(building_relay).vertices(v,2),buildings(building_relay).vertices(next_v,1),buildings(building_relay).vertices(next_v,2));
                    if intersection==1
                        %There is an intersection
                        map_indoor_distances(i,j)=min(norm(posBS-point),map_indoor_distances(i,j));
                    end
                end
           end        
        end
    end
end


