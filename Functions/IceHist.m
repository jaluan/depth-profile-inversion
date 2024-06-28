%% This function returns the distributions of ice histories that are 
% accepted, and within 1/2 sd of measured concentrations

function [uval,uval2,uval3,bval,bval2,bval3]=IceHist(model,dO_i,bt) 

uval = []; uval2 = []; uval3 = [];
for nw = 1:model.Nwalk
    I = find(model.walker{nw}.status == 1); %accepted
    I2 = find(model.walker{nw}.status == 1 & model.walker{nw}.restot<=4); %2sd
    I3 = find(model.walker{nw}.status == 1 & model.walker{nw}.restot<=1); %1sd
    if (~isempty(I))
        uval = [uval(:);model.walker{nw}.up(1,I)'];
    end
    if (~isempty(I2))
        uval2 = [uval2(:);model.walker{nw}.up(1,I2)'];
    end
    if (~isempty(I3))
        uval3 = [uval3(:);model.walker{nw}.up(1,I3)'];
    end
end

% Burial in % since 1 Ma
bval=interp1(dO_i,bt,uval); bval2=interp1(dO_i,bt,uval2); bval3=interp1(dO_i,bt,uval3);
