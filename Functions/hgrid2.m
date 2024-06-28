%% This function calculates grid for exhumation plot for second sample
function histgrid=hgrid2(tbin,zint,tint,model)  

histgrid = zeros(size(tbin));
dzi = zint(2)-zint(1);

n0 = model.Mmp+model.Nsmp; %changed to second sample
z0 = model.z0;
Nt = length(tint);
Nz = length(zint);

for nw=1:model.Nwalk
    I = find(model.walker{nw}.status == 1);
    for i=1:length(I)

        T1 = model.walker{nw}.u(2,I(i));
        z1 = model.walker{nw}.u(n0+1,I(i));
        dT2 = model.walker{nw}.u(n0+2,I(i));
        dz2 = 10^model.walker{nw}.u(n0+3,I(i));
        dT3 = model.walker{nw}.u(n0+4,I(i));
        dz3 = 10^model.walker{nw}.u(n0+5,I(i));
        E4 = 10^model.walker{nw}.u(n0+6,I(i));

        T2 = T1 + dT2;
        z2 = z1 + dz2;
        T3 = T2 + dT3;
        z3 = z2 + dz3;
        T4 = model.age; %this requires age > Tdg+dT2+dT3
        z4 = z3 + (model.age - T3)*E4;

        Tm = [0,T1,T2,T3,T4];
        zm = [0,z1,z2,z3,z4];

        zinterp = interp1(Tm,zm,tint,'linear',2*z0);        
        izs = 1+floor(zinterp/dzi); %bin index at each tsfine
        for itsfine=1:Nt
            if izs(itsfine)<=Nz
                histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
            end
        end
    end        
end
histgrid = histgrid/max(histgrid(:));