function compile_synthetic_data(scn)

close all; %clear
ns = 1; %number of samples

for i=1:ns %loop over samples

    % get CN properties
    CNprop = getCNprop;

    %% Define 'Sample' details 
    sample{i}.name = ['SyntSample' num2str(i)];
    sample{i}.type = 'bedrock';
    sample{i}.batchid = 'none';
    sample{i}.Nnuc = 2;
    sample{i}.lat=59.8482; %latitude (dd)
    sample{i}.lon=8.6600; %longitude (dd)
    sample{i}.elevation = 1714; %elevation (m)

    % Calculate atm pressure (hPa) from elevation
    p = ERA40atm(sample{i}.lat,sample{i}.lon,sample{i}.elevation);    
    sample{i}.pressure = p; 

    sample{i}.minDgla = 10e-3; %minimum deglaciation age (Ma)
    sample{i}.maxDgla = 20e-3; %maximum deglaciation age (Ma)
    sample{i}.shield=1; %topographic shielding factor [0-1]
    sample{i}.thick=2; %sample thickness (cm), used to correct surface production below
    sample{i}.density = CNprop.rho;%density in g/cm3
    rho = sample{i}.density; %Density of rock sample (g/cm3), for local use
    
    %% Production %%
    % Spallation attenuation and thickness correction
    % Lspal=attenuationlength(sample{i}.lat,sample{i}.lon,sample{i}.elevation,p); %Calculated from CronusCalc functions based on site cutoff rigidity, not considering terrain shielding %[g/cm2]
    Lspal = 155; %Alternative, constant value [g/cm2]
    sample{i}.production.Lspal=Lspal; %[g/cm2]
    
    %Thickness correction, spallation
    sf_spal = exp(-sample{i}.thick/2*rho/Lspal); %Factor to correct production 
    % rate for thickness of sample, sets surface production = production midway 
    % through sample. Make sure this is not already factored in to site-specific 
    % production rates. Set to 1 to exclude.

    %Define depths below surface z/rho cm/(g/cm3) [g/cm^2] for fitting of production profiles
    D_m = 100; %Depth, changed below (rho)
    z_m = linspace(0,10,100);
    z_D = D_m*z_m.^3/10*rho; %denser depth-grid near surface
    maxZ = 5300; %maxdepth (g/cm2) for muon-production profile used for fitting
    % of exponentials below. If this depth is very large, exponential terms will
    % be dominated by fast muon production, which isn't ideal. 1200g/cm2=4.5m
    % with rho ~2.65-2.7, 5300 g/cm2=20 m, Test effect of this choice

     %% 10Be, mandatory
    % Spallation surface production in atoms / g qtz / year
    % 4.01 = Be10-spallation at SLHL from Borchers et al.,
    % 2016, Table 7, for the st scaling framework
    sample{i}.production.P10spal = 4.01.*stone2000(sample{i}.lat,p,1); 
    sample{i}.production.P10spal = sample{i}.production.P10spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction

    %Muon production following Balco 2017, 
    %Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
    mc10.k_neg = 0.00191 .* 0.704 .* 0.1828; %summary cross-section for negative muon capture (at/muon)
    mc10.sigma0 = 0.280e-30; %x-section for fast muon production at 1 Gev
    mc10.Natoms = 2.006e22; %Oxygen atoms pr gram Quartz
    % Fit muon production profile calculated with P_mu_total_alpha1.m 
    % with three exponential terms following Balco 2017, Eq. 7/Fig. 15
    p10_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,3,mc10,0);
    sample{i}.production.P10_Lm1=p10_muons.L(1); %attenuation, first exponential term
    sample{i}.production.P10_Lm2=p10_muons.L(2); %attenuation, second exponential term
    sample{i}.production.P10_Lm3=p10_muons.L(3); %attenuation, third exponential term
    sample{i}.production.P10_m1 = p10_muons.P(1); %production, first exponential term
    sample{i}.production.P10_m2 = p10_muons.P(2); %production, second exponential term
    sample{i}.production.P10_m3 = p10_muons.P(3); %production, third exponential term
    shield_fac10_m1 = exp(-sample{i}.thick/2*rho/p10_muons.L(1)); %sample thickness correction
    shield_fac10_m2 = exp(-sample{i}.thick/2*rho/p10_muons.L(2)); %sample thickness correction
    shield_fac10_m3 = exp(-sample{i}.thick/2*rho/p10_muons.L(3)); %sample thickness correction
    sample{i}.production.P10_m1 = sample{i}.production.P10_m1*shield_fac10_m1*sample{i}.shield; % production, first exponential term, corrected for sample thickness and topographic shielding
    sample{i}.production.P10_m2 = sample{i}.production.P10_m2*shield_fac10_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
    sample{i}.production.P10_m3 = sample{i}.production.P10_m3*shield_fac10_m3*sample{i}.shield; %production, third exponential term, corrected for sample thickness and topographic shielding
  
    %Muon production following Heisinger 2002
    sample{i}.production.P10_Lfm = CNprop.Lfm; %fast muons
    sample{i}.production.P10_Lnmc = CNprop.Lnmc; %negative muons
    P10tot= 4.01.*stone2000(sample{i}.lat,p,1)/(1-CNprop.pr_fm_Be-CNprop.pr_nmc_Be); 
    sample{i}.production.P10_fm = CNprop.pr_fm_Be*P10tot; %fast muons
    sample{i}.production.P10_nmc = CNprop.pr_nmc_Be*P10tot; %negative muon capture
    shield_fac10_fm = exp(-sample{i}.thick/2*rho/CNprop.Lfm); %sample thickness correction
    shield_fac10_nmc = exp(-sample{i}.thick/2*rho/CNprop.Lnmc); %sample thickness correction
    sample{i}.production.P10_fm = sample{i}.production.P10_fm*shield_fac10_fm*sample{i}.shield; % production, first exponential term, corrected for sample thickness and topographic shielding
    sample{i}.production.P10_nmc = sample{i}.production.P10_nmc*shield_fac10_nmc*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
 
    %% 26Al, mandatory
    % Spallation surface production in atoms / g qtz / year
    % 27.93 = Al26-spallation at SLHL from Borchers et al.,
    % 2016, Table 7, for the st scaling framework
    sample{i}.production.P26spal = 27.93.*stone2000(sample{i}.lat,p,1);
    sample{i}.production.P26spal = sample{i}.production.P26spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
    
    %Muon production following Balco 2017, 
    %Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
    mc26.k_neg = 0.0133 .* 0.296 .* 0.6559; %summary cross-section for negative muon capture (at/muon)
    mc26.sigma0 = 3.89e-30; %x-section for fast muon production at 1 Gev
    mc26.Natoms = 1.003e22; %Si atoms pr gram Quartz
    % Fit muon production profile calculated with P_mu_total_alpha1.m 
    % with three exponential terms following Balco 2017, Eq. 7/Fig. 15
    p26_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,3,mc26,0);
    sample{i}.production.P26_Lm1=p26_muons.L(1); %attenuation, first exponential term
    sample{i}.production.P26_Lm2=p26_muons.L(2); %attenuation, second exponential term
    sample{i}.production.P26_Lm3=p26_muons.L(3); %attenuation, third exponential term
    sample{i}.production.P26_m1 = p26_muons.P(1); %production, first exponential term
    sample{i}.production.P26_m2 = p26_muons.P(2); %production, second exponential term
    sample{i}.production.P26_m3 = p26_muons.P(3); %production, third exponential term
    shield_fac26_m1 = exp(-sample{i}.thick/2*rho/p26_muons.L(1)); %sample thickness correction
    shield_fac26_m2 = exp(-sample{i}.thick/2*rho/p26_muons.L(2)); %sample thickness correction
    shield_fac26_m3 = exp(-sample{i}.thick/2*rho/p26_muons.L(3)); %sample thickness correction
    sample{i}.production.P26_m1 = sample{i}.production.P26_m1*shield_fac26_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding
    sample{i}.production.P26_m2 = sample{i}.production.P26_m2*shield_fac26_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
    sample{i}.production.P26_m3 = sample{i}.production.P26_m3*shield_fac26_m3*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
   
    %Muon production following Heisinger 2002
    sample{i}.production.P26_Lfm = CNprop.Lfm; %fast muons
    sample{i}.production.P26_Lnmc = CNprop.Lnmc; %negative muons
    P26tot= 27.93.*stone2000(sample{i}.lat,p,1)/(1-CNprop.pr_fm_Al-CNprop.pr_nmc_Al); 
    sample{i}.production.P26_fm = CNprop.pr_fm_Al*P26tot; %fast muons
    sample{i}.production.P26_nmc = CNprop.pr_nmc_Al*P26tot; %negative muon capture
    shield_fac26_fm = exp(-sample{i}.thick/2*rho/CNprop.Lfm); %sample thickness correction
    shield_fac26_nmc = exp(-sample{i}.thick/2*rho/CNprop.Lnmc); %sample thickness correction
    sample{i}.production.P26_m1 = sample{i}.production.P26_fm*shield_fac26_fm*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding
    sample{i}.production.P26_m2 = sample{i}.production.P26_nmc*shield_fac26_nmc*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
    
%     %% 21Ne, optional, needs updating
%     if sample{i}.Nnuc == 3            
%         % 21Ne/10Be production ratio from Balco and Shuster, 2009:
%         % 4.08 +- 0.37 at /g qtz /yr, multiplied w. 10Be spal. prod
%         sample{i}.production.P21spal = 4.08*4.01.*stone2000(sample{i}.lat,p,1);
%         sample{i}.production.P21spal = sample{i}.production.P21spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
% 
%         % Muon production following Balco 2017, 
%         % Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
%         % Muons 21Ne. These are not calibrated, taken from Fernandez-Mosquera 2010
%         mc21.Natoms = 1.0003e22; %Si atoms pr gram Quartz
%         mc21.k_neg = 0.296.*0.6559.*0.0029; %summary cross-section for negative muon capture (at/muon)
%         mc21.sigma190 = 0.79e-27; %x-section for fast muon production at 1 Gev
% 
%         % Fit muon production profile calculated with P_mu_total_alpha1.m 
%         % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
%         p21_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc21,0);
%         sample{i}.production.P21_Lm1=p21_muons.L(1); %attenuation, first exponential term
%         sample{i}.production.P21_Lm2=p21_muons.L(2); %attenuation, second exponential term
%         sample{i}.production.P21_m1 = p21_muons.P(1); %production, first exponential term
%         sample{i}.production.P21_m2 = p21_muons.P(2); %production, second exponential term
%         shield_fac21_m1 = exp(-sample{i}.thick/2*rho/p21_muons.L(1)); %sample thickness correction
%         shield_fac21_m2 = exp(-sample{i}.thick/2*rho/p21_muons.L(2)); %sample thickness correction
%         sample{i}.production.P21_m1 = sample{i}.production.P21_m1*shield_fac21_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding 
%         sample{i}.production.P21_m2 = sample{i}.production.P21_m2*shield_fac21_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding 
%     end

    %% Set model parameters %%
%     d18Op = 3.9;  % d18O threshold (3.5-5)
%     Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
%     Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
%     dT2   = 0.26;  % Time change from Tdgla [Myr] (0-2 Myr)
%     dZ2   = -0.8;   % Depth change during T2 [log10(m)] (-1 to 1)
%     dT3   = 0.14;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
%     dZ3   = 0.5;  % Depth change during T3 [log10(m)] (-1 to 1)
%     E0    = 0.99;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)
   [d18Op,Tdgla,Z1,dT2,dZ2,dT3,dZ3,E0] = scenarios(scn); %load parameter values

    up = [d18Op,Tdgla,Z1,dT2,dZ2,dT3,dZ3,E0]; %Pack model parameters to vector structure
%     zp=[0.0,Z1,Z1+10^(dZ2),Z1+10^(dZ2)+10^(dZ3),Z1+10^(dZ2)+10^(dZ3)+(3-(Tdgla+dT2+dT3))*10^(E0)];
%     tp=[0.0,Tdgla*10^3,(Tdgla+dT2)*10^3,(Tdgla+dT2+dT3)*10^3,3000];
%     model_fig(d18Op,zp,tp)
    model.synpmval=up;

    model.Nsnr = 1; %number of samples
    model.Nfree = 2; %Number of free model depth points
    model.data{i} = sample{i};
    model.Nsmp = 2*model.Nfree + 2; %number of sample specific parameters
    model.data{i}.depths=(0:0.4:2)';%[0;0.5;1.5]; %depths of datapoints (m)
    model.Ndp = length(model.data{1}.depths); %Number of data points in depth profile
    model.Nnc = 2; %number of nuclides
    model.Nds = model.Nnc*model.Ndp; %number of data per sample (nuclides*depths)
    model.age = 3.0; %Max time (Myr)
    model.z0 = 20; %Max depth
    model.Temp = 1.0;
    model.Mmp = 2; %number of generic model parameters
    model.data{i}.density=CNprop.rho;
%     model.data.production=sample{i}.production;

    %load and prepare d18O data
    load('d18Ocurves.mat','Age','d18O_4ky');
    I = find(Age < model.age*1e6);
    model.dO = d18O_4ky(I);
    model.dOt = Age(I);

    %% Run forward model to get 'true' nuclide concentrations
    [gm] = forward_bedrockvJ4(up,model,CNprop); %Balco 2007 muons
    N10uncT = getCosmoUnc(gm(1:model.Ndp),10,'Be10',0); %calculates likely error in % on 10Be
    N10uncT=max(N10uncT,2.5*ones(size(N10uncT))); %Set minimum error to 2.5%
    err10true = gm(1:model.Ndp).*N10uncT/100; %Corresponding error in at/g
    N26uncT = getCosmoUnc(gm(model.Ndp+1:end),10,'Al26',0); %calculates likely error in % on 26Al
    N26uncT=max(N26uncT,2.5*ones(size(N26uncT))); %Set minimum error to 2.5%
    err26true = gm(model.Ndp+1:end).*N26uncT/100; %Corresponding error in at/g

    % Draw a random number around 'true' concentration with 2.5% std dev,
    % NB: not used in manuscript.
    for j=1:100
        gmr(j,:)=gm+0.025*gm.*randn(size(gm));
    end

    % Calculate uncertainty of measurements based on 'historical data' in Aarhus lab + AMS
    N10unc = getCosmoUnc(gmr(:,1:model.Ndp),10,'Be10',0); %calculates likely error in % on 10Be
    N10unc=max(N10unc,2.5*ones(size(N10unc))); %Set minimum error to 2.5%
    err10 = gmr(:,1:model.Ndp).*N10unc/100; %N10spal.*N10unc/100; Corresponding error in at/g on synthetic samples
    N26unc = getCosmoUnc(gmr(:,model.Ndp+1:end),10,'Al26',0); %calculates likely error in % on 26Al
    N26unc=max(N26unc,2.5*ones(size(N26unc))); %Set minimum error to 2.5%
    err26 = gmr(:,model.Ndp+1:end).*N26unc/100; %N26spal.*N26unc/100; Corresponding error in at/g on synthetic samples

%     figure(), hold on
%     tiledlayout(1,3)
%     set(gcf,'units','pixels','position',[1500,1500,1000,600]);
%     nexttile
%     patch([gmr(1,1:model.Ndp)-err10(1,1:model.Ndp) fliplr(gmr(1,1:model.Ndp)+err10(1,1:model.Ndp))],...
%         [model.data{1}.depths' flipud(model.data{1}.depths)'],0.5*[1 1 1],'FaceAlpha',0.3), hold on
%     plot(gmr(1,1:model.Ndp),model.data{1}.depths','ko')
%     plot(gm(1:model.Ndp),model.data{1}.depths','rx')
%     xlabel('N10 (at/g)'), ylabel('Depth (m)')
%     set(gca,'ydir','reverse','Fontsize',20)
%     legend('Uncertainty','Random','Balco true','Location','southeast')
%     nexttile
%     patch([gmr(1,model.Ndp+1:end)-err26(1,1:model.Ndp) fliplr(gmr(1,model.Ndp+1:end)+err26(1,1:model.Ndp))],...
%         [model.data{1}.depths' flipud(model.data{1}.depths)'],0.5*[1 1 1],'FaceAlpha',0.3), hold on
%     plot(gmr(1,model.Ndp+1:end),model.data{1}.depths','ko')
%     plot(gm(1,model.Ndp+1:end),model.data{1}.depths','rx')
%     xlabel('N26 (at/g)'), ylabel('Depth (m)')
%     set(gca,'ydir','reverse','Fontsize',20)

    sample{i}.depths=model.data{1}.depths;
    sample{i}.N10 = gmr(:,1:model.Ndp); %N10spal(1);
    sample{i}.dN10 = err10;
    sample{i}.N10true = gm(1:model.Ndp); %N10spal(1);
    sample{i}.dN10true = err10true;
    sample{i}.N26 = gmr(:,model.Ndp+1:end); %N26spal(1);
    sample{i}.dN26 = err26;
    sample{i}.N26true = gm(model.Ndp+1:end); %N26spal(1);
    sample{i}.dN26true = err26true;
    sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
    sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);

%     if sample{i}.Nnuc == 3, % this section needs updating if 21Ne is included
%         sample{i}.N21=num(ij:ij+sample{i}.Ndp-1,17);
%         sample{i}.dN21=num(ij:ij+sample{i}.Ndp-1,18);
%         
%         %compute ratios
%         sample{i}.r1021 = sample{i}.N10./sample{i}.N21;
%         sample{i}.dr1021 = sample{i}.r1021.*sqrt((sample{i}.dN21./sample{i}.N21).^2+(sample{i}.dN10./sample{i}.N10).^2);
%         sample{i}.r2621 = sample{i}.N26./sample{i}.N21;
%         sample{i}.dr2621 = sample{i}.r2621.*sqrt((sample{i}.dN21./sample{i}.N21).^2+(sample{i}.dN10./sample{i}.N10).^2);
%     end        
      
    nexttile
    patch([sample{i}.r2610(1,:)-sample{i}.dr2610(1,:) fliplr(sample{i}.r2610(1,:)+sample{i}.dr2610(1,:))],...
        [model.data{1}.depths' flipud(model.data{1}.depths)'],0.5*[1 1 1],'FaceAlpha',0.3), hold on
    plot(sample{i}.r2610(1,:),model.data{1}.depths','ko')
    % plot(gm(model.Ndp+1:end)./gm(1:model.Ndp),model.data{1}.depths','rx')
    xlabel('N26/N10 (at/g)'), ylabel('Depth (m)')
    set(gca,'ydir','reverse','Fontsize',20)

end
% savefile = ['data/synthetic_data_fw_dp_v4_scn', num2str(scn),'_n.mat'];
% save(savefile,'model','sample');