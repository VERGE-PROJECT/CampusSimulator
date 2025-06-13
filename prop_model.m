function L = prop_model(d,model,f,hUT,hBS,hE,NLOS_condition,indoor_condition,d2Dindoor)
%PROP_MODEL Summary of this function goes here
%   Detailed explanation goes here
%d: distance in m
%model: 0:UMi, 1:UMa, 2: D2D_WINNER, 3: D2D_Xia
%f: frequency in GHz
%hUT,hBS: Height of UE and BS in m
%NLOS_condition: 0:LOS, 1:NLOS
%indoor_condition: 0: outdoor, 1: indoor
%d2Dindoor: distance inside the building (in m)



switch model
    case 0  %UMi model
        d3D=sqrt(d*d+(hBS-hUT)*(hBS-hUT));
        dBP=4*(hBS-hE)*(hUT-hE)*f*1E9/3E8;
        %First compute path loss under LOS conditions:
        if d<=dBP
            L_LOS=32.4+21*log10(d3D)+20*log10(f);
        else
            L_LOS=32.4+40*log10(d3D)+20*log10(f)-9.5*log10(dBP*dBP+(hBS-hUT)*(hBS-hUT));
        end
        
        if NLOS_condition==0
            %LOS
            L=L_LOS;
        else
            %NLOS
            L_NLOS=35.3*log10(d3D)+22.4+21.3*log10(f)-0.3*(hUT-1.5);
            L=max(L_NLOS,L_LOS);
        end
        if indoor_condition==1
           %Add penetration losses: 
           Lindoor=5-10*log10(0.3*power(10,-0.1*(2+0.2*f))+0.7*power(10,-0.1*(5+4*f)))+0.5*d2Dindoor; %Low-loss model
           %Lindoor=5-10*log10(0.7*power(10,-0.1*(23+0.3*f))+0.3*power(10,-0.1*(5+4*f)))+0.5*d2Dindoor; %High-loss model
           L=L+Lindoor;
        end
       

    case 4  %D2D model (Siemens model 3GPP based on pedestrian test environment as described in UMTS30.03 and in 3GPP TR25.942
        if NLOS_condition==0
            %LOS model: Free-space loss
            L=32.44+20*log10(f)+20*log10(max(d,1)); %Assume that if distance is less than 1m the path loss at 1m is considered.
        else
            %NLOS model
            L=19+30*log10(f)+40*log10(max(d,1));
        end
        if indoor_condition==1
           %Add penetration losses: 
           Lindoor=5-10*log10(0.3*power(10,-0.1*(2+0.2*f))+0.7*power(10,-0.1*(5+4*f)))+0.5*d2Dindoor; %Low-loss model
           %Lindoor=5-10*log10(0.7*power(10,-0.1*(23+0.3*f))+0.3*power(10,-0.1*(5+4*f)))+0.5*d2Dindoor; %High-loss model
           L=L+Lindoor;
        end
      
end


end

