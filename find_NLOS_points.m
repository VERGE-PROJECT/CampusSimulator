function map_NLOS = find_NLOS_points(sizeX,sizeY,pixel_size,posBS,buildings,map_indoor)

%find_NLOS_points: This function generates a map of dimensions
%sizeX,sizeY in pixels of size "pixel_size". The map includes the NLOS
%points for the BS in position BS considering the buildings.
%   Detailed explanation goes here
%- sizeX,sizeY: size of the scenario in pixels
%- pixel_size;
%- posBS: Position of the BS (or the transmitter in case of D2D)
%%- buldings: each component buildings(n) is a structure where vertices
%includes a matrix with the points of all vertices defining the building in
%m
%- map_indoor: map with the positions that correspond to indoor points


map_NLOS=zeros(sizeX,sizeY);
num_buildings=size(buildings,2);
for i=1:sizeX
    for j=1:sizeY
        pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
        
        %Method: 
        %If it is an outdoor point -> it is NLOS if there is at least 1
        %intersection with the building segments.
        %If it is an indoor point -> it is NLOS if there are at least 2
        %intersections with the building segments.
        
        if map_indoor(i,j)==1
            max_intersections=2;
        else
            max_intersections=1;
        end
        
        num_intersect=0;
        for b=1:num_buildings
           num_vertices=size(buildings(b).vertices,1);
           for v=1:num_vertices
               next_v=mod(v,num_vertices)+1;
               [aux,~]=intersect_lines(posBS(1),posBS(2),pos(1),pos(2),buildings(b).vertices(v,1),buildings(b).vertices(v,2),buildings(b).vertices(next_v,1),buildings(b).vertices(next_v,2));
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

