% Extract modelled depth at time TE1 (Ma) and modelled times at depth ZT1 (m)
function [E1,TZ1] = interpModel(model,TE1,ZT1)

n0=2;
E1 = []; TZ1 = [];

for nw = 1:model.Nwalk %loop over walkers
    I = find(model.walker{nw}.status > 0); %accepted walkers
    if (~isempty(I))
        
        E1w = [];
        TZ1w = [];
        
        for k=1:length(I)
                    
            T1 = model.walker{nw}.u(2,I(k));
            z1 = model.walker{nw}.u(n0+1,I(k));
            dT2 = model.walker{nw}.u(n0+2,I(k));
            dz2 = 10^model.walker{nw}.u(n0+3,I(k));
            dT3 = model.walker{nw}.u(n0+4,I(k));
            dz3 = 10^model.walker{nw}.u(n0+5,I(k));
            E4 = 10^model.walker{nw}.u(n0+6,I(k));
        
            T2 = T1 + dT2;
            z2 = z1 + dz2;
            T3 = T2 + dT3;
            z3 = z2 + dz3;
            T4 = model.age; %this requires age > Tdg+dT2+dT3
            z4 = z3 + (model.age - T3)*E4;
            
            Tm = [0,T1,T2,T3,T4]; %Ma
            zm = [0,z1,z2,z3,z4]; %m

            E1w(k) = interp1(Tm,zm,TE1); %interpolate time-depth vector at time TE1
            
            if (z4 > ZT1)
                TZ1w(k) = interp1(zm,Tm,ZT1); %interpolate depth-time vectors at depth ZT1
            else
                TZ1w(k) = (ZT1 - z4)/E4 + T4; %Calculate time if ZT1 is below max sample depth
            end
           
        end
        
        %save things
        E1 = [E1(:);E1w(:)];
        TZ1 = [TZ1(:);TZ1w(:)];
        
    end
end
