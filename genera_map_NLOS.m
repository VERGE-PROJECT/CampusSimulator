function map_NLOS = genera_map_NLOS(sizeX,sizeY,pixel_size,resolution_dist_NLOS,dist_max_NLOS,posBS,model,hUT,f,num_aleat)

%GENERA_MAP_NLOS Summary of this function goes here
%   Detailed explanation goes here
%sizeX,sizeY: size of the scenario in pixels
%pixel_size;
%resolution_dist_NLOS: Resolution of the NLOS distances.
%dist_max_NLOS: Maximum distance for NLOS.
%posBS: Position of the BS (or the transmitter in case of D2D)
%model: Propagation model (for computing LOS probability)
%hUT: Height of the UE
%f: Frequency (GHz). Only used for computing LOS probability in D2D.
%num_aleat: Integer number to randomize the execution of this function (leave 0 to have it deterministic).


distances_NLOS=[resolution_dist_NLOS:resolution_dist_NLOS:dist_max_NLOS]';
angles_NLOS=zeros(size(distances_NLOS,1),360); %360: angle with granularity 1º
d_index=1;
d_index_max=size(distances_NLOS,1);
%start_angle=1;
for d=distances_NLOS'
    
    pNLOS=1-LOS_prob(d,model,hUT,f);
    range=round(pNLOS*360);
    
    %p_NLOS_values(d_index,1)=pNLOS;
    %range_values(d_index,1)=range;
    
    %The following process marks a total of "range" angles in angles_NLOS.
    %This is done in groups of consecutive "size_block" angles starting
    %from "start_angle". After a group the next start position is shifted
    %110º further. In case a group cannot be completed because some
    %positions are already marked from previous distances, it continues
    %marking in the first available position.
    
    aux=sum(angles_NLOS(d_index,:)); %These are the angles already marked as NLOS coming from previous distances
    start_angle=mod((d_index-1)*40+num_aleat,360)+1;
    size_block=20;
    
    while aux<range
        angle=start_angle;
        num_pos_completed=0;
        endloop=0;
        while endloop==0
            if angles_NLOS(d_index,angle)==0
                angles_NLOS(d_index:d_index_max,angle)=1; %NLOS at distance d_index will also be NLOS at larger distances
                num_pos_completed=num_pos_completed+1;
                aux=aux+1;
                if (num_pos_completed==size_block) || (aux==range)
                    endloop=1; 
                end
            end
            angle=mod(angle,360)+1;
            if (angle==start_angle)
                endloop=1; %There are not enough positions to fill. 
                aux=range; %To get out of the other loop
            end
        end
        start_angle=mod(angle+110+num_aleat,360)+1;  %Start at another position.
    end
    d_index=d_index+1;
end

%Build map of NLOS positions
map_NLOS=ones(sizeX,sizeY);

pixels_max=ceil(dist_max_NLOS/pixel_size);
posBS_index=ceil(posBS/pixel_size);


for i=max(1,posBS_index(1)-pixels_max):min(sizeX,posBS_index(1)+pixels_max)
    for j=max(1,posBS_index(2)-pixels_max):min(sizeY,posBS_index(2)+pixels_max)
        pos=pixel_size*[i-0.5,j-0.5];  %Position is the center of the pixel
        dist=norm(pos-posBS);
        if dist<=dist_max_NLOS
            dist_index=floor(dist/resolution_dist_NLOS);
            angle_index=1+mod(round(atan2(pos(2)-posBS(2),pos(1)-posBS(1))*(180/pi)+360),360); %angle of atan2 is between -pi and pi, we round it and convert it to be between 1 and 360
        
            if dist_index==0
                map_NLOS(i,j)=0;
            else
                map_NLOS(i,j)=angles_NLOS(dist_index,angle_index);
            end

        end
    end
end

