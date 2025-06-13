function prob = LOS_prob(d,model,hUT,f)
%LOS_PROB Summary of this function goes here
%   Detailed explanation goes here
%d: distance in m
%model: 0: UMi, 1:UMa
%hUT: Height of terminal in m. %Only used by UMa
%f: Frequency (GHz)
switch model
    case 0  %UMi model
        if d<=18
            prob=1;
        else
            prob=(18/d)+(1-18/d)*exp(-d/36);
        end
        
    case 1  %UMa model
        if d<=18
            prob=1;
        else
            if hUT<=13
                C=0;
            else
                C=power((hUT-13)/10,1.5);
            end
            
            prob=((18/d)+(1-18/d)*exp(-d/63))*(1+C*(5/4)*power(d/100,3)*exp(-d/150));
        end

        
    case 4  %D2D model (Siemens model 3GPP based on pedestrian test environment as described in UMTS30.03 and in 3GPP TR25.942)
        xmin=power(10,(13.44-10*log10(f))/20);
        xmax=50;
        if d<xmin
            prob=1;
        elseif d<xmax
            prob=1-(d-xmin)/(xmax-xmin);
        else
            prob=0;
        end

end

end

