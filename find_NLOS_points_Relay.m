function map_NLOS = find_NLOS_points_Relay(sizeX,sizeY,pixel_size,posBS,buildings,map_buildings,dmax)

%find_NLOS_points_Ralay: This function generates a map of dimensions
%sizeX,sizeY in pixels of size "pixel_size". The map includes the NLOS
%points for the relay in a given position considering the buildings.
%   Detailed explanation goes here
%- sizeX,sizeY: size of the scenario in pixels
%- pixel_size;
%- posBS: Position of the Relay
%%- buldings: each component buildings(n) is a structure where vertices
%includes a matrix with the points of all vertices defining the building in
%m
%- map_buildings: map with the building of each point
%- dmax: maximum distance for NLOS.(points at d>dmax will always be NLOS)


map_NLOS=ones(sizeX,sizeY);  %By default all points are NLOS
num_buildings=size(buildings,2);

pixels_max=ceil(dmax/pixel_size);
posBS_index=ceil(posBS/pixel_size);

building_relay=map_buildings(posBS_index(1),posBS_index(2));
%Filter only the relevant buildings to check:
if building_relay>0
    %Relay indoor: Only relevant building is the building where it is
    %located.
    num_valid_buildings=1;
    buildings_aux(num_valid_buildings)=buildings(building_relay);   
else
   %Relay outdoor: Only relevant buildings are those with at least one vertex in distance closer than dmax.
   num_valid_buildings=0;
   for b=1:num_buildings
        num_vertices=size(buildings(b).vertices,1);
        for v=1:num_vertices
            dist=norm(buildings(b).vertices(v,:)-posBS);
            if dist<=dmax
                %valid building:
                num_valid_buildings=num_valid_buildings+1;
                buildings_aux(num_valid_buildings)=buildings(b);
                %Stop the loop of vertices
                break
            end
        end
   end
end

%Check only a square centered in the position of the relay and distance
%2*pixels_max. All the other points are NLOS.
for i=max(1,posBS_index(1)-pixels_max):min(sizeX,posBS_index(1)+pixels_max)
    for j=max(1,posBS_index(2)-pixels_max):min(sizeY,posBS_index(2)+pixels_max)
        pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
        

        %Method: 

        %If distance between point and relay is >d_max -> point is NLOS in any case.
        %If relay is indoor:
        %  -> all outdoor points are NLOS.
        %  -> an indoor point of the same building is NLOS if there is at
        %  least 1 intersection with the same building. (Note that this
        %  will only occur in case of very irregular buildings).
        
        %If relay is outdoor
        %  -> all indoor points are NLOS.
        %  -> an outdoor point is NLOS if there is at least 1 intersection
        %  with the valid buildings
        
        %Based on this, we only need to check the points of the same
        %building as the relay (considering building=0 if points are
        %outdoor), and then check the intersections with the relevant
        %buildings.
        
        %Building where the point (i,j) is located (0 if outdoor).
        building_point=map_buildings(i,j);
        
        if building_relay==building_point
            %Both the relay and the point are at the same building or both of them are outdoor (if building_relay==0):
            
            if norm(pos-posBS)<=dmax
            
                %Check intersections with the set of valid buildings:
                num_intersect=0;
                max_intersections=1;
            
                map_NLOS(i,j)=0;  %Initially assume LOS and check if it has to be NLOS.
                for b=1:num_valid_buildings
                    num_vertices=size(buildings_aux(b).vertices,1);
                    for v=1:num_vertices
                        next_v=mod(v,num_vertices)+1;
                        [aux,~]=intersect_lines(posBS(1),posBS(2),pos(1),pos(2),buildings_aux(b).vertices(v,1),buildings_aux(b).vertices(v,2),buildings_aux(b).vertices(next_v,1),buildings_aux(b).vertices(next_v,2));
                        num_intersect=num_intersect+aux;
                        if num_intersect>=max_intersections
                            %NLOS point
                            map_NLOS(i,j)=1;
                            %end the loop of v
                            break 
                        end
                    end
                    if map_NLOS(i,j)==1
                        break
                    end
                end
            end
        end
    end
end

