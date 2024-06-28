
function [N10Out,N26Out] = cosmoStepDepth(N10In,N26In,dt,pf,prod,erate,zd,ss)

%************* Cosmo parameters **************
P10spal = prod.P10spal; %atoms/(g*yr)
P10m1 = prod.P10_m1;
P10m2 = prod.P10_m2;
P10m3 = prod.P10_m3;

P26spal = prod.P26spal; %atoms/(g*yr)
P26m1 = prod.P26_m1;
P26m2 = prod.P26_m2;
P26m3 = prod.P26_m3;

Lspal = prod.Lspal; % Attenuation lengths g/cm^2
P10Lm1 = prod.P10_Lm1;
P10Lm2 = prod.P10_Lm2;
P10Lm3 = prod.P10_Lm3;
P26Lm1 = prod.P26_Lm1;
P26Lm2 = prod.P26_Lm2;
P26Lm3 = prod.P26_Lm3;

TBe = 1.387e6;
TAl = 0.705e6;
lambda_Be = log(2)/TBe;
lambda_Al = log(2)/TAl;
rho = 2.65; %density g/cm3


%********************************************
switch ss
    case 'ss' %steady state erosion, full exposure
        %spallation
        fBe = lambda_Be + rho*erate*100/Lspal;
        fAl = lambda_Al + rho*erate*100/Lspal;
        N10 = P10spal*exp(-rho*100.*zd./Lspal)/fBe; %N10 = P10spal/fBe; 
        N26 = P26spal*exp(-rho*100.*zd./Lspal)/fAl; %N26 = P26spal/fAl;
        %1st muon pathway
        fBe = lambda_Be + rho*erate*100/P10Lm1;
        fAl = lambda_Al + rho*erate*100/P26Lm1;
        N10 = N10 + P10m1*exp(-rho*100.*zd./P10Lm1)/fBe; %N10 = N10 + P10m1/fBe;
        N26 = N26 + P26m1*exp(-rho*100.*zd./P26Lm1)/fAl; %N26 = N26 + P26m1/fAl;
        %2nd muon pathway
        fBe = lambda_Be + rho*erate*100/P10Lm2;
        fAl = lambda_Al + rho*erate*100/P26Lm2;
        N10 = N10 + P10m2*exp(-rho*100.*zd./P10Lm2)/fBe; %N10 = N10 + P10m2/fBe;
        N26 = N26 + P26m2*exp(-rho*100.*zd./P26Lm2)/fAl; %N26 = N26 + P26m2/fAl;
        %3rd muon pathway
        fBe = lambda_Be + rho*erate*100/P10Lm3;
        fAl = lambda_Al + rho*erate*100/P26Lm3;
        N10 = N10 + P10m3*exp(-rho*100.*zd./P10Lm3)/fBe; %N10 = N10 + P10m3/fBe;
        N26 = N26 + P26m3*exp(-rho*100.*zd./P26Lm3)/fAl; %N26 = N26 + P26m3/fAl;


    case 'eb' %concentration change given erate, dt, exposure

        P10z_spal = pf*P10spal*exp(-rho*100*zd/Lspal);
        P10z_m1 = pf*P10m1*exp(-rho*100*zd/P10Lm1);
        P10z_m2 = pf*P10m2*exp(-rho*100*zd/P10Lm2);
        P10z_m3 = pf*P10m3*exp(-rho*100*zd/P10Lm3);
        P26z_spal = pf*P26spal*exp(-rho*100*zd/Lspal);
        P26z_m1 = pf*P26m1*exp(-rho*100*zd/P26Lm1);
        P26z_m2 = pf*P26m2*exp(-rho*100*zd/P26Lm2);
        P26z_m3 = pf*P26m3*exp(-rho*100*zd/P26Lm3);
        
        N10 = N10In*exp(-lambda_Be*dt); % decay since last timestep
        fBe = lambda_Be+rho*erate*100/Lspal;
        N10 = N10 + P10z_spal'*(1.0-exp(-fBe*dt))/fBe;
        fBe = lambda_Be+rho*erate*100/P10Lm1;
        N10 = N10 + P10z_m1'*(1.0-exp(-fBe*dt))/fBe;
        fBe = lambda_Be+rho*erate*100/P10Lm2;
        N10 = N10 + P10z_m2'*(1.0-exp(-fBe*dt))/fBe;
        fBe = lambda_Be+rho*erate*100/P10Lm3;
        N10 = N10 + P10z_m3'*(1.0-exp(-fBe*dt))/fBe;
        
        N26 = N26In*exp(-dt*lambda_Al); % decay since last timestep
        fAl = lambda_Al+rho*erate*100/Lspal;
        N26 = N26 + P26z_spal'*(1.0-exp(-fAl*dt))/fAl;
        fAl = lambda_Al+rho*erate*100/P26Lm1;
        N26 = N26 + P26z_m1'*(1.0-exp(-fAl*dt))/fAl;
        fAl = lambda_Al+rho*erate*100/P26Lm2;
        N26 = N26 + P26z_m2'*(1.0-exp(-fAl*dt))/fAl;
        fAl = lambda_Al+rho*erate*100/P26Lm3;
        N26 = N26 + P26z_m3'*(1.0-exp(-fAl*dt))/fAl;
            
end

%Output
N10Out = N10;
N26Out = N26;