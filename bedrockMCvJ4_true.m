function bedrockMCvJ4_true(snr,fnr,rnr)
% vJ2 updated Feb. 2020 by JLA includes the possibility of extra data points at depth
% vJ2_1: testing Balco 2017 muon implementation using
% forward_bedrock_vJ2_1.m
% vJ3: Aug 2021, better documentation, and stone-calibrated production rates from
% compile_data_vJ3
% vJ4: Nov 2022, three exponentials to muon fit to depth of app. 20 m
% vJ4_true: inverts 'true' values of synthetic depth profiles (no randomn
% noise)

% Inverse MCMC model exploring complex exposure and exhumation histories of 
% samples with measured cosmogenic nuclides from (formerly) glaciated
% landscapes. Delineates the ensemble of parameter values with lowest 
% data misfit.

% Inputs: 
% snr: scalar or vector containing sample numbers of samples in input data
% file generated from compile_mn_data.m. If multiple numbers are given, the
% exposure history of the samples will be assumed to be the same, while
% exhumation may vary. Note that the function retrieves last deglaciation 
% timing from the first sample only (should be similar).
% fnr: filenumber, controls which prescribed scenario to run
% rnr: random error, not in use.

% Written by David Lundbek Egholm
% Modified by Jane Lund Andersen to also include samples at depth and 
% (optional) 21Ne measurements.

close all;

%% Cosmogenic  halflives
CNprop = getCNprop;

%% load data
if fnr==21
    load('./data/Lamuvaara_v4_surf_n.mat','sample'), dpname = 'Lamuvaara_v4_surf_n'; % 'sample', 'model'
elseif fnr==22
    load('./data/Naakakarhakka_v4_surf_n.mat','sample'), dpname = 'Naakakarhakka_v4_surf_n'; % 'sample', 'model'
elseif fnr==23
    load('./data/Gausta_v4_surf.mat','sample'), dpname = 'Gausta_v4_surf'; % 'sample', 'model'
elseif fnr==24
    load('./data/Lysefjord_v4_surf.mat','sample'), dpname = 'Lysefjord_v4_surf'; % 'sample', 'model'
elseif fnr==25
    load('./data/Karmøy_v4_surf.mat','sample'), dpname = 'Karmøy_v4_surf'; % 'sample', 'model'
elseif fnr==26
    load('./data/Andøya_v4_surf.mat','sample'), dpname = 'Andøya_v4_surf'; % 'sample', 'model'
else
    load(['./data/synthetic_data_fw_dp_v4_scn',num2str(fnr),'_n'],'sample','model')
    dpname = ['Syn_v4_scn',num2str(fnr),'_true_n'];
    % load(['./data/synthetic_data_fw_dp_v4_scn',num2str(fnr),'a_n'],'sample','model') %only surface sample
    % dpname = ['Syn_v4_scn',num2str(fnr),'_truesurf_n'];
end

%% set model parameters
%number of samples
model.Nsnr = length(snr);

%Number of free depth points
model.Nfree = 2;

%number of sample specific parameters
model.Nsmp = 2*model.Nfree + 2;

%initialize
% models = struct();

model.age = 3.0; %Max time (Myr)
model.z0 = 20; %Max depth
model.Temp = 1.0;
model.Mmp = 2; %number of generic model parameters
model.mp{1}.name = ['d18O threshold'];
model.mp{1}.vmin = 3.5;
model.mp{1}.vmax = 5.0;
model.mp{2}.name = ['Time of deglaciation (Tdg)'];
model.mp{2}.vmin = sample{1}.minDgla; %10e-3; %Myr
model.mp{2}.vmax = sample{1}.maxDgla; %20e-3;

%load and prepare d18O data
load('d18Ocurves.mat','Age','d18O_4ky');
I = find(Age < model.age*1e6);
model.dO = d18O_4ky(I);
model.dOt = Age(I);

%loop number of samples and for model parameters
for i=1:model.Nsnr
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.name =  ['Z at Tdg, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmax = 0.05; %m
    model.mp{model.Mmp+(i-1)*model.Nsmp+2}.name =  ['dT2, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmax = 1.5; %Myr
    model.mp{model.Mmp+(i-1)*model.Nsmp+3}.name =  ['dz2, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmin = -1; %0.1 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmax = 1; %log 10 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+4}.name =  ['dT3, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmax = 1.5; %Myr 
    model.mp{model.Mmp+(i-1)*model.Nsmp+5}.name =  ['dz3, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmin = -1;%0.1 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmax = 1; %log 10 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+6}.name =  ['E4, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmin = -2; %log10 -2 =0.01m/Myr
    model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmax = 2; %2=100 m/Myr

    %add sample info to models
    model.data{i} = sample{snr(i)};
    
    %depths of datapoints (m)
    model.data{i}.depths=sample{i}.depths;

    %Number of data points in depth profile
    model.Ndp = length(sample{i}.depths); %Nb: for now only handling multiple depth profiles of same length and #nuclides

    %number of nuclides
    model.Nnc = sample{i}.Nnuc;

    %total number of data per sample (nuclides*depths)
    model.Nds = model.Nnc*model.Ndp;

end

%save some general MCMC paramers
model.Nwalk = 10; %6; Number of walkers
model.burnin = 10e3; %5e3; Length of burnin
model.Nmod = 100e3; %50e3; Number of models
model.Nmax = 1e6; %500000; Maximum number of models

%number of total model parameters
model.Nmp = model.Mmp + model.Nsnr*model.Nsmp;

for i=1:model.Nmp
    umin(i) = model.mp{i}.vmin;
    umax(i) = model.mp{i}.vmax;
    du0(i) = model.mp{i}.vmax - model.mp{i}.vmin;
end
umin = umin(:);
umax = umax(:);

%data and covariance
dobs=zeros(1,model.Nds*model.Nsnr); sigd=zeros(1,model.Nds*model.Nsnr); %initiate data observation and error vectors
if fnr > 20 %Real data
    for i=1:model.Nsnr
        for j=1:model.Ndp
            dobs((i-1)*model.Nds+j) = model.data{i}.N10(j);
            dobs((i-1)*model.Nds+j+model.Ndp) = model.data{i}.N26(j);
            sigd((i-1)*model.Nds+j) = model.data{i}.dN10(j);
            sigd((i-1)*model.Nds+j+model.Ndp) = model.data{i}.dN26(j);
            if model.Nnc == 3
               dobs((i-1)*model.Nds+j+2*model.Ndp) = model.data{i}.N21(j);
               sigd((i-1)*model.Nds+j+2*model.Ndp) = model.data{i}.dN21(j);
            end
        end
    end
else %synthetic data
    for i=1:model.Nsnr
        for j=1:model.Ndp
            dobs((i-1)*model.Nds+j) = model.data{i}.N10true(j);
            dobs((i-1)*model.Nds+j+model.Ndp) = model.data{i}.N26true(j);
            sigd((i-1)*model.Nds+j) = model.data{i}.dN10true(j);
            sigd((i-1)*model.Nds+j+model.Ndp) = model.data{i}.dN26true(j);
            if model.Nnc == 3
               dobs((i-1)*model.Nds+j+2*model.Ndp) = model.data{i}.N21true(j);
               sigd((i-1)*model.Nds+j+2*model.Ndp) = model.data{i}.dN21true(j);
            end
        end
    end
end

Cobs = model.Temp*diag(sigd.^2);
Cobsinv = inv(Cobs);

%%%%%% Initiate random models generation %%%%%%%

%Initialize random sequence
rng('default');

%walker starting points
if model.Nwalk > 1
    for i=1:model.Nmp
        wini(i,:) = 0.8*(randperm(model.Nwalk)-1)/(model.Nwalk-1)+0.1;
    end
else
    wini = randi([0,4],model.Nmp,model.Nwalk)/4;
end

%loop walkers
for nw = 1:model.Nwalk

    %walker starting point - for initial parameter vector
    
    for i = 1:model.Nmp
        u(i) = (1-wini(i,nw))*model.mp{i}.vmin + wini(i,nw)*model.mp{i}.vmax;
    end
                
 
    %initialize
    minres = 1e20; %not used
    res_current = 1e20; %current residual - first model run
    restot = 0;
    acount = 0;
    bcount = 0;
    rcount = 0;
    accrat = 0;
    status = zeros(model.Nmax,1);
    erosion_rec = zeros(model.Nmax,1);
    up_rec = zeros(model.Nmp,model.Nmax);
    u_rec = zeros(model.Nmp,model.Nmax);
    N10_rec = zeros(model.Nmax,1);
    N26_rec = zeros(model.Nmax,1);
    restot_rec = zeros(model.Nmax,1);
    accrat_rec = zeros(model.Nmax,1);
    k_rec = zeros(model.Nmax,1);
    duR = zeros(model.Nmp,1);

   
    accfac = 1e-2;
    ktarget = 0.025; %not used
    
    %run models
    mi = 0; %model iteration
    k = 0.01; %initial step length

    while ((mi < model.Nmax)&&(acount < model.Nmod))
        
        mi = mi + 1;

%         disp(['nw = ',num2str(nw),'/',num2str(model.Nwalk),' mi = ',num2str(mi),'/',num2str(model.Nmax),' bcount = ',num2str(bcount),' acount = ',num2str(acount),' accrat = ',num2str(accrat),' k = ',num2str(k)]);
    
       %***** step length ******
       
       
       %acceptance ratio
        if (mi > 100)
            accrat = (sum(abs(status((mi-100):(mi-1))))+1)/100;
        elseif (bcount < model.burnin)
            accrat = 0.1;
        else 
            accrat = 0.3;
        end
        
        %burnin
        if (bcount < 0.5*model.burnin) 
            
            k = k*((1-accfac) + accfac*accrat/0.1);
           
            if ((mi > 100)&&(bcount < 2)) k = 1.0;
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi;
            end
                       
            
            
        elseif (bcount < model.burnin) 

            k = k*((1-accfac) + accfac*accrat/0.2);
           
            if ((mi > 100)&&(bcount < 2)) k = 1.0;
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi;
            end
            
            %k = 2*ktarget;
        
        elseif (acount < model.Nmod)

            k = k*((1-accfac) + accfac*accrat/0.3);    
            
            %k = k*(accrat + 0.1)/(0.3 + 0.1)
            
            %k = ktarget;
            
        end   
        
        if (k > 0.5) k = 0.5; end
        
        if (bcount < model.burnin) model.Temp = 1.0 + 10.0*(model.burnin-bcount)/model.burnin;
        else model.Temp = 1.0;
        end
        
        %********* propose new parameters ************
        
         %random step
        du = 0.5*randn(model.Nmp,1).*du0(:);
        
        %proposed model
        up = u(:) + k*du(:);
        
        %retake
        while ((any(up(:) < umin(:)))||(any(up(:) > umax(:)))||((up(2)+up(4)+up(6)) > model.age)) %note last clause may need updating if model parameters change
        
            %random step
            du = 0.5*randn(model.Nmp,1).*du0(:);

            %proposed model
            up = u(:) + k*du(:);
            
        end
        
        
        %********** Forward model *****************
        if model.Nnc(1) == 2 %10Be, 26Al 
            [gmp] = forward_bedrockvJ4(up,model,CNprop);   
        elseif model.Nnc(1) == 3 %21Ne
            [gmp] = forward_bedrockvJ3_Ne(up,model,CNprop);  %update
        end
    
               
        %Acceptance criteria
        %restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:));
        restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:))/model.Temp;
        rfrac = exp(-0.5*restot)/exp(-0.5*res_current);
        alpha = rand(1);
    
    
        %if model is accepted
        if ((alpha < rfrac)||(mi == 1))
        
            u = up;
            gm = gmp;
            res_current = restot;
        
            %accepted after burnin
            if (bcount > model.burnin)
        
                status(mi) = 1;
                acount = acount + 1;
                
            %if accepted during burnin    
            else
        
                status(mi) = -1;
                bcount = bcount + 1;
            
            end
        
        %rejected
        else
        
            status(mi) = 0;
            
            %rejected
            rcount = rcount + 1;                 
    
        end
       
        %save things
        up_rec(:,mi) = up(:);
        u_rec(:,mi) = u(:);
        gm_rec(:,mi) = gm(:);
        restot_rec(mi) = res_current;
        accrat_rec(mi) = accrat;
        k_rec(mi) = k;
                
    end

    %change status flag for models that are rejected during burnin
    Imin=find(status == 1,1);  % find index of first accepted model after burnin
    I = find(status(1:Imin) == 0); %find indices of rejected models in burnin
    status(I) = -2; %set status of these models to -2
    
    model.walker{nw}.status = status(1:mi);
    model.walker{nw}.up = up_rec(:,1:mi);
    model.walker{nw}.u = u_rec(:,1:mi);
    model.walker{nw}.gm = gm_rec(:,1:mi);
    model.walker{nw}.restot = restot_rec(1:mi);
    model.walker{nw}.acount = acount;
    model.walker{nw}.bcount = bcount; 
    model.walker{nw}.rcount = rcount;
    model.walker{nw}.accrate = accrat_rec(1:mi);
    model.walker{nw}.kstep = k_rec(1:mi);

end

%save output

str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

savefile = ['models/' dpname '.mat'];
save(savefile,'model');