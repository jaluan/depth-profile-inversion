% Extract true depth at time TE1 (Ma) and true time at depth ZT1 (m)
function [E1t,TZ1t] = trueModel(synpmval,modelAge,TE1,ZT1)

n0=2;

T1 = synpmval(1,2);
z1 = synpmval(1,n0+1);
dT2 = synpmval(1,n0+2);
dz2 = 10^synpmval(1,n0+3);
dT3 = synpmval(1,n0+4);
dz3 = 10^synpmval(1,n0+5);
E4 = 10^synpmval(1,n0+6);

T2 = T1 + dT2;
z2 = z1 + dz2;
T3 = T2 + dT3;
z3 = z2 + dz3;
T4 = modelAge; %this requires age > Tdg+dT2+dT3
z4 = z3 + (T4 - T3)*E4;

Tm = [0,T1,T2,T3,T4]; %Ma
zm = [0,z1,z2,z3,z4]; %m

E1t = interp1(Tm,zm,TE1); %interpolate time-depth vector at time TE1

if (z4 > ZT1)
    TZ1t = interp1(zm,Tm,ZT1); %interpolate depth-time vectors at depth ZT1
else
    TZ1t = (ZT1 - z4)/E4 + T4; %Calculate time if ZT1 is below max sample depth
end

end

