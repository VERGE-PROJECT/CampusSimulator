function map_indoor_distances= find_indoor_distances(sizeX,sizeY,pixel_size,posBS,buildings,map_buildings)
%FIND_INDOOR_DISTANCES Summary of this function goes here
%   This function computes the 2D indoor distance between a BS and a point
%   that is inside a building.
%Inputs:
%- Scenario parameters (sizeX,sizeY) in pixels, pixel_size in m
%- posBS: position of the BS
%- buildings: struct with the vertices of each building
%- map_buildings: Map that indicates the building where each pixel belongs
%(0 if it is an outdoor point)

map_indoor_distances=zeros(sizeX,sizeY);
max_distance_in_scenario=sqrt((sizeX*pixel_size)^2+(sizeY*pixel_size)^2);
for i=1:sizeX
    for j=1:sizeY
        pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
        
        b=map_buildings(i,j);
        if b>0
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
           map_indoor_distances(i,j)=max_distance_in_scenario;
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


