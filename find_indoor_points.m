function [map_indoor_points,map_buildings] = find_indoor_points(sizeX,sizeY,pixel_size,buildings)
%FIND_INDOOR_POINTS Summary of this function goes here
% This function takes as input a set of buildings characterizes by their
% vertices in a scenario of size X, sizeY and pixel_size. Then, it provides
% a map including the positions that correspond to indoor points. It also
% provides a map with the positions of each building.
%   Detailed explanation goes here
%sizeX: number of pixels in horizontal dimension
%sizeY: number of pixels in vertical dimension
%pixel_size: m
%buldings: each component buildings(n) is a structure where vertices
%includes a matrix with the points of all vertices defining the building in
%m


num_buildings=size(buildings,2);
map_indoor_points=zeros(sizeX,sizeY);
map_buildings=zeros(sizeX,sizeY);

for i=1:sizeX
    for j=1:sizeY
        %if i==326 && j==140
        %    fprintf("Ara.\n");
        %end
        pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
        pos=pos+0.01*[rand,rand]; %It is convenient to introduce a certain randomness to avoid that the intersections fall exactly in the vertices
        

        for b=1:num_buildings
           num_intersect=0;
           num_vertices=size(buildings(b).vertices,1);
           for v=1:num_vertices
              next_v=mod(v,num_vertices)+1;  %for the last vertex the segment is between this vertex and the first one.
              [aux,~]=intersect_lines(0,0,pos(1),pos(2),buildings(b).vertices(v,1),buildings(b).vertices(v,2),buildings(b).vertices(next_v,1),buildings(b).vertices(next_v,2));
              num_intersect=num_intersect+aux;
           end
           
           if mod(num_intersect,2)==1
               %Point is inside building b. 
               map_indoor_points(i,j)=1;
               map_buildings(i,j)=b;
               %Stop the loop.
               break
           end  
        end
    end
end



end

