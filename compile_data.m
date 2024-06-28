function compile_data()
% 
% vJ2 updated Feb. 2020 by JLA: includes the possibility of extra data points at depth
% vJ3 updated Aug. 2021 by JLA: stone-calibrated production rates, Balco
% 2017 muons and Hesinger muons, 21Ne optional
% vJ4: Use three exponentials to fit depth profiles to 20 m depth
% 
% This code reads sample data from an excelfile, calculates production
% parameters using 'St' scaling, and saves a data structure for use in
% bedrockMC codes.
% 
% Note: at present this and subsequent codes are set up to handle 10Be and
% 26Al (mandatory)
% 
% Input excelfile format:
% Each sample is presented with a line for each sample in depth profile. 
% The first line contains sample site information (Sample name,
% lat (dd), lon (dd), elevation (m), topographic shielding factor (0-1), 
% sample thickness (cm) and density (g/cm3), year of sampling (yr), timing 
% of last deglaciation and uncertainty (ka), number of nuclides measured, 
% number of depth points in profile (determines how many lines to skip 
% before ereading next sample. 
% Each line contains the sample depth below surface (cm), measured nuclide
% concentrations and uncertainties in this order: N10, dN10, N26, dN26,
% N21(optional), dN21(optional)

% Cited literature:
% Balco 2017: Production rate calculations for cosmic-ray-muon-produced 
% 10Be and 26Al benchmarked against geological calibration data. Quaternary
% Geochronology 39 (150-173).
% Borchers et al., 2016: Borchers, B., Marrero, S., Balco, G., Caffee, M., 
% Goehring, B., Lifton, N., Nishiizumi, K., Phillips, F., Schaefer, J. and 
% Stone, J., 2016. Geological calibration of spallation production rates in 
% the CRONUS-Earth project. Quaternary Geochronology, 31, pp.188-198.
% Fernandez-Mosquera 2010:
% Heisinger 2002: 
% Balco & Shuster 2009: Balco, G. and Shuster, D.L., 2009. Production rate
% of cosmogenic 21Ne in quartz estimated from 10Be, 26Al, and 21Ne 
% concentrations in slowly eroding Antarctic bedrock surfaces. Earth and 
% Planetary Science Letters, 281(1-2), pp.48-58.

close all; clear
addpath('Functions')

 % Get CN properties (half-lives etc.)
CNprop = getCNprop;

%% Read data info from Excelfile
ns = 2; %number of samples/sites, multiple depthpoints from one core profile equals 1 sample. 
% Set ns>1 for joined inversion of multiple depth profiles with similar exposure history.

% Loop over samples
for i=1:ns
%% Read data from excel file, this needs modification if depth profiles are inverted separately.
if i==1
    excelfile = 'data/Naakakarhakka.xlsx'; %dpname = 'Naakakarhakka_v4';
elseif i==2
    excelfile = 'data/Lamuvaara.xlsx'; %dpname = 'Lamuvaara_v4';
end
dpname = 'Naa_Lam_v4';

[num,text,~] = xlsread(excelfile);
    

% %% Only surface sample
% excelfile = 'data/Naakakarhakka.xlsx'; dpname = 'Naakakarhakka_v4_surf';
% excelfile = 'data/Lamuvaara.xlsx'; dpname = 'Lamuvaara_v4_surf';
% [num,text,~] = xlsread(excelfile,"Sheet2");

%% 
    ij=1; %row start

    sample{i}.Ndp = num(ij,11); %Number of data points at depth for sample i
    sample{i}.Nnuc = num(ij,10); %number of nuclides
    
    sample{i}.name = text(ij+1,1); %Sample ID
    
    sample{i}.lat=num(ij,1); %latitude (dd)
    sample{i}.lon=num(ij,2); %longitude (dd)
    sample{i}.elevation = num(ij,3); %elevation (m)
    
    % Calculate atm pressure (hPa) from elevation
    p = ERA40atm(sample{i}.lat,sample{i}.lon,sample{i}.elevation);    
    sample{i}.pressure = p; 
    
    sample{i}.shield=num(ij,4); %topographic shielding factor [0-1]
    sample{i}.thick=num(ij,5); %sample thickness (cm), used to correct surface production below
    sample{i}.density = num(ij,6); %density in g/cm3
    rho = sample{i}.density; %Density of rock sample (g/cm3), for local use
    sample{i}.sampleyr = num(ij,7); %sample year, not in use
    
    % timing of last deglaciation
    sample{i}.minDgla = (num(ij,8)-num(ij,9))*1e-3; %last deglaciation minimum estimate, [Myr]
    sample{i}.maxDgla = (num(ij,8)+num(ij,9))*1e-3; %last deglaciation maximum estimate, [Myr]
  
    %Read depth specific data
    sample{i}.depths = num(ij:ij+sample{i}.Ndp-1,12)/100; %sample depths converted to meters
    sample{i}.N10 = num(ij:ij+sample{i}.Ndp-1,13); %N10 measured
%     sample{i}.dN10 = num(ij:ij+sample{i}.Ndp-1,14); %err10
    sample{i}.dN10 = num(ij:ij+sample{i}.Ndp-1,19); %err10 % min 2.5 %
    sample{i}.N26 = num(ij:ij+sample{i}.Ndp-1,15); %N26 measured
%     sample{i}.dN26 = num(ij:ij+sample{i}.Ndp-1,16); %err26
    sample{i}.dN26 = num(ij:ij+sample{i}.Ndp-1,20); %err26 % min 2.5 %
    
    %calculate N26/N10 ratio
    sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
    sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
    
    %Site-specific production parameters

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
    % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
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
    % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
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
    sample{i}.production.P26_m3 = sample{i}.production.P26_m3*shield_fac26_m3*sample{i}.shield; %production, third exponential term, corrected for sample thickness and topographic shielding

    %% 21Ne, optional
    if sample{i}.Nnuc == 3
        sample{i}.N21=num(ij:ij+sample{i}.Ndp-1,17);
        sample{i}.dN21=num(ij:ij+sample{i}.Ndp-1,18);
        
        %compute ratios
        sample{i}.r1021 = sample{i}.N10./sample{i}.N21;
        sample{i}.dr1021 = sample{i}.r1021.*sqrt((sample{i}.dN21./sample{i}.N21).^2+(sample{i}.dN10./sample{i}.N10).^2);
        sample{i}.r2621 = sample{i}.N26./sample{i}.N21;
        sample{i}.dr2621 = sample{i}.r2621.*sqrt((sample{i}.dN21./sample{i}.N21).^2+(sample{i}.dN10./sample{i}.N10).^2);
            
        % 21Ne/10Be production ratio from Balco and Shuster, 2009:
        % 4.08 +- 0.37 at /g qtz /yr, multiplied w. 10Be spal. prod
        sample{i}.production.P21spal = 4.08*4.01.*stone2000(sample{i}.lat,p,1);
        sample{i}.production.P21spal = sample{i}.production.P21spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction

        % Muon production following Balco 2017, 
        % Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
        % Muons 21Ne. These are not calibrated, taken from Fernandez-Mosquera 2010
        mc21.Natoms = 1.0003e22; %Si atoms pr gram Quartz
        mc21.k_neg = 0.296.*0.6559.*0.0029; %summary cross-section for negative muon capture (at/muon)
        mc21.sigma190 = 0.79e-27; %x-section for fast muon production at 1 Gev

        % Fit muon production profile calculated with P_mu_total_alpha1.m 
        % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
        p21_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,3,mc21,0);
        sample{i}.production.P21_Lm1=p21_muons.L(1); %attenuation, first exponential term
        sample{i}.production.P21_Lm2=p21_muons.L(2); %attenuation, second exponential term
        sample{i}.production.P21_Lm3=p21_muons.L(3); %attenuation, third exponential term
        sample{i}.production.P21_m1 = p21_muons.P(1); %production, first exponential term
        sample{i}.production.P21_m2 = p21_muons.P(2); %production, second exponential term
        sample{i}.production.P21_m3 = p21_muons.P(3); %production, third exponential term
        shield_fac21_m1 = exp(-sample{i}.thick/2*rho/p21_muons.L(1)); %sample thickness correction
        shield_fac21_m2 = exp(-sample{i}.thick/2*rho/p21_muons.L(2)); %sample thickness correction
        shield_fac21_m3 = exp(-sample{i}.thick/2*rho/p21_muons.L(3)); %sample thickness correction
        sample{i}.production.P21_m1 = sample{i}.production.P21_m1*shield_fac21_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding 
        sample{i}.production.P21_m2 = sample{i}.production.P21_m2*shield_fac21_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding 
        sample{i}.production.P21_m3 = sample{i}.production.P21_m3*shield_fac21_m3*sample{i}.shield; %production, third exponential term, corrected for sample thickness and topographic shielding
    end
    %% Look for next sample in this row
    ij=ij+sample{i}.Ndp; 
end

savefile = ['data/' dpname '_n.mat'];
save(savefile,'sample','excelfile');