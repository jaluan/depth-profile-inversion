function [gm] = forward_bedrockvJ4(up,model,CNprop)

% function [gm] = forward_bedrockvJ4(up,model,CNprop)
% Forward model for transient CN integration of bedrock samples 

% Inputs: up = model parameter vector, model = model setup (from bedrockMC), 
% CNprop = structure with halflives and decay constants of radioactive 
% nuclides (from getCNprop.m)

% Output: gm = predicted data vector

% Written by David Lundbek Egholm
% Modified by Jane Lund Andersen to also include samples at depth and
% updated production handling

%% model parameters

%generic glaciation parameters
d18Op = up(1);
T1 = up(2);

%load and prepare d18O data
%% load and prepare d18O data
dO = model.dO;
dOt = model.dOt;

%time step
model.dt = 1000; %yr

%loop samples
for i=1:model.Nsnr

    n0 = (i-1)*model.Nsmp+model.Mmp; %parameter number start
    
    z1 = up(n0+1);
    dT2 = up(n0+2); %Myr
    dz2 = 10^up(n0+3); %m
    dT3 = up(n0+4);
    dz3 = 10^up(n0+5); %m
    E4 = 10^up(n0+6); %m/Myr

    T2 = T1 + dT2;
    z2 = z1 + dz2;
    
    T3 = T2 + dT3;
    z3 = z2 + dz3;
    
    T4 = model.age; %this requires age > Tdg+dT2+dT3
    z4 = z3 + (model.age - T3)*E4;
    
    mT = [0,T1,T2,T3,T4]*1e6; %Myr to yr
    mz = [0,z1,z2,z3,z4];
    
    % Check starting condition of model
    if (z4 > model.z0) % Depth at start of model greater than max depth
        maxtime = interp1(mz,mT,model.z0);
        ssbc = 0;
    else
        maxtime = model.age*1e6;
        ssbc = 1;
    end
    
    nt = ceil(maxtime/model.dt); 
    time = maxtime*linspace(0,1,nt);
    
    
    %**************************************
    % Exhumation history
    %**************************************
%     burial = interp1(mT,mz,time); %Depths at times in 'time' %old
    
    Ndp=model.Ndp; %Number of data points in depth profile
    depths=model.data{i}.depths; %depths of datapoints
    
    burial = zeros(Ndp,nt); %initiate 'burial' matrix with depths at times in 'time', for all data sample points in depth profile
    burial(1,:) = interp1(mT,mz,time); %Surface depths at times in 'time'
    
    if Ndp > 1 % if there are samples at depth
        burial(2:end,:) = burial(1,:)+depths(2:end).*ones(Ndp-1,nt); %add their depths below surface to 'burial'
    end
    
    dOn = interp1(dOt,dO,time,'linear',dO(end)); %d18O values at 'time'
    pfac = ones(nt,1); %controls exposure, modified below
    pfac(dOn > d18Op) = 0; %No exposure when d18O values above threshold (d18Op)
    pfac(time < 25e3) = 0; %correct exposure around deglaciation: set no
                           %exposure last 25 kyr
    pfac(time < T1*1e6) = 1; %then add exposure from time of deglaciation
                                
    %CN production parameters
    rho = model.data{i}.density;
    Lspal = model.data{i}.production.Lspal;
    P10spal = model.data{i}.production.P10spal;
    P26spal = model.data{i}.production.P26spal;

    % Be-10
    P10m1 = model.data{i}.production.P10_m1;
    P10Lm1 = model.data{i}.production.P10_Lm1;
    P10m2 = model.data{i}.production.P10_m2;
    P10Lm2 = model.data{i}.production.P10_Lm2;
    P10m3 = model.data{i}.production.P10_m3;
    P10Lm3 = model.data{i}.production.P10_Lm3;
    % Al-26
    P26m1 = model.data{i}.production.P26_m1;
    P26Lm1 = model.data{i}.production.P26_Lm1;
    P26m2 = model.data{i}.production.P26_m2;
    P26Lm2 = model.data{i}.production.P26_Lm2;
    P26m3 = model.data{i}.production.P26_m3;
    P26Lm3 = model.data{i}.production.P26_Lm3;

    
    N10 = zeros(Ndp,nt); % changed '1' to 'Ndp'
    N26 = zeros(Ndp,nt); % --
    
    if (ssbc == 0) %depth of sample at model start greater than maxdepth
        
        N10(:,nt) = 0.0; % changed '(nt)' to '(:,nt)'
        N26(:,nt) = 0.0; % --

    else
        %assume steady state concentration at starting point
        
        erate = E4*1e-6;
        
        %spallation
        fBe = CNprop.lambda_Be + rho*erate*100/Lspal;
        fAl = CNprop.lambda_Al + rho*erate*100/Lspal;
        N10(:,nt) = P10spal*exp(-rho*100*burial(:,nt)/Lspal)/fBe; % changed '(nt)' to '(:,nt)'
        N26(:,nt) = P26spal*exp(-rho*100*burial(:,nt)/Lspal)/fAl; % changed '(nt)' to '(:,nt)'
        %1st muon pathway
        fBe = CNprop.lambda_Be + rho*erate*100/P10Lm1;
        fAl = CNprop.lambda_Al + rho*erate*100/P26Lm1;
        N10(:,nt) = N10(:,nt) + P10m1*exp(-rho*100*burial(:,nt)/P10Lm1)/fBe; % changed '(nt)' to '(:,nt)'
        N26(:,nt) = N26(:,nt) + P26m1*exp(-rho*100*burial(:,nt)/P26Lm1)/fAl; % changed '(nt)' to '(:,nt)'
        %2nd muon pathway
        fBe = CNprop.lambda_Be + rho*erate*100/P10Lm2;
        fAl = CNprop.lambda_Al + rho*erate*100/P26Lm2;
        N10(:,nt) = N10(:,nt) + P10m2*exp(-rho*100*burial(:,nt)/P10Lm2)/fBe; % changed '(nt)' to '(nt,:)'
        N26(:,nt) = N26(:,nt) + P26m2*exp(-rho*100*burial(:,nt)/P26Lm2)/fAl; % changed '(nt)' to '(nt,:)'
        %3rd muon pathway
        fBe = CNprop.lambda_Be + rho*erate*100/P10Lm3;
        fAl = CNprop.lambda_Al + rho*erate*100/P26Lm3;
        N10(:,nt) = N10(:,nt) + P10m3*exp(-rho*100*burial(:,nt)/P10Lm3)/fBe; % changed '(nt)' to '(nt,:)'
        N26(:,nt) = N26(:,nt) + P26m3*exp(-rho*100*burial(:,nt)/P26Lm3)/fAl; % changed '(nt)' to '(nt,:)'

    end
    
        
    %integrate time
    for kk=(nt-1):-1:1
        
        dt = (time(kk+1)-time(kk));
        pf = .5*(pfac(kk+1)+pfac(kk));
        bz = burial(:,kk); % changed '(kk)' to '(:,kk)', Depths at time step kk
        erate = 100*(burial(1,kk+1)-burial(1,kk))/dt; % changed '(kk)' to '(1,kk)'
        
        P10z_spal = pf*P10spal*exp(-rho*100*bz/Lspal);
        P10z_m1 = pf*P10m1*exp(-rho*100*bz/P10Lm1);
        P10z_m2 = pf*P10m2*exp(-rho*100*bz/P10Lm2);
        P10z_m3 = pf*P10m3*exp(-rho*100*bz/P10Lm3);
        P26z_spal = pf*P26spal*exp(-rho*100*bz/Lspal);
        P26z_m1 = pf*P26m1*exp(-rho*100*bz/P26Lm1);
        P26z_m2 = pf*P26m2*exp(-rho*100*bz/P26Lm2);
        P26z_m3 = pf*P26m3*exp(-rho*100*bz/P26Lm3);
              
        N10(:,kk) = N10(:,kk+1)*exp(-dt*CNprop.lambda_Be); % changed '(kk)' to '(:,kk+1)'
        ff = CNprop.lambda_Be+rho*erate/Lspal;
        N10(:,kk) = N10(:,kk) + P10z_spal*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        ff = CNprop.lambda_Be+rho*erate/P10Lm1;
        N10(:,kk) = N10(:,kk) + P10z_m1*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        ff = CNprop.lambda_Be+rho*erate/P10Lm2;
        N10(:,kk) = N10(:,kk) + P10z_m2*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        ff = CNprop.lambda_Be+rho*erate/P10Lm3;
        N10(:,kk) = N10(:,kk) + P10z_m3*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
        N26(:,kk) = N26(:,kk+1)*exp(-dt*CNprop.lambda_Al); % changed '(kk)' to '(:,kk+1)'
        ff = CNprop.lambda_Al+rho*erate/Lspal;
        N26(:,kk) = N26(:,kk) + P26z_spal*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        ff = CNprop.lambda_Al+rho*erate/P26Lm1;
        N26(:,kk) = N26(:,kk) + P26z_m1*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        ff = CNprop.lambda_Al+rho*erate/P26Lm2;
        N26(:,kk) = N26(:,kk) + P26z_m2*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        ff = CNprop.lambda_Al+rho*erate/P26Lm3;
        N26(:,kk) = N26(:,kk) + P26z_m3*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
    end

    gm=zeros(1,model.Nds*model.Nsnr);
    for j=1:Ndp
        gm((i-1)*model.Nds+j) = N10(j,1);
        gm((i-1)*model.Nds+Ndp+j) = N26(j,1);
    end
end

