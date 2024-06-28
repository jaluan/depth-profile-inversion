% This script produces Fig. 1 - illustrating the inherent burial-erosion
% trade off for 10Be-26Al surface samples

close all; clear
set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
% addpath('Data')
% addpath('Functions')

%% Get production rates and sample characteristics of synthetic sample
load data/synthetic_data_fw_dp_v4_scn1_n.mat sample
CNprop = getCNprop;

%% Define synthetic scenarios
Nscn=2; %number of scenarios

scn=1; %Scenario A
d18Op(scn) = 4.65; %4.55 % d18O threshold (3.5-5)
Tdgla(scn) = 5e-3; % Time of last deglaciation [Myr] (2-15e-3)
Z1(scn)    = 0.0001; %0.0185;    % Depth at time of deglaciation [m] (0-0.05)
dT2(scn)   = 0.1;  % Time change from Tdgla [Myr] (0-2 Myr)
dZ2(scn)   = -0.05; %0.17;   % Depth change during T2 [log10(m)] (-1 to 1)
dT3(scn)   = 1.9;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
dZ3(scn)   = -0.5;  % Depth change during T3 [log10(m)] (-1 to 1)
E0(scn)    = -1.2; %0.1;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

scn=2; %Scenario B
d18Op(scn) = 3.75; %3.55 % d18O threshold (3.5-5)
Tdgla(scn) = 11e-3; % Time of last deglaciation [Myr] (2-15e-3)
Z1(scn)    = 0.0001; %0.0185; % Depth at time of deglaciation [m] (0-0.05)
dT2(scn)   = 0.7;  % Time change from Tdgla [Myr] (0-2 Myr)
dZ2(scn)   = -0.95; %-0.8;   % Depth change during T2 [log10(m)] (-1 to 1)
dT3(scn)   = 1.28;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
dZ3(scn)   = -0.95; %-0.8;  % Depth change during T3 [log10(m)] (-1 to 1)
E0(scn)    = 0.42;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

%% colours and line styles
mcols=[.3 .5 .7; .5 .4 .3; .7 .9 .5; .4 .2 .3; .5 .7 .2]; 
mstyls={'o','s','v','^','d','x'};
lcols=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840]; %Matlab default
lstyls={'-','-.',':','--','-'};

%% Choose model parameters
model.Nsnr = 1; %number of samples
model.Nfree = 2; %Number of free model depth points
model.data{1} = sample{1};
model.Nsmp = 2*model.Nfree + 2; %number of sample specific parameters
model.data{1}.depths=(0:0.2:2)';%[0;0.5;1.5]; %depths of datapoints (m)
model.Ndp = length(model.data{1}.depths); %Number of data points in depth profile
model.Nnc = 2; %number of nuclides
model.Nds = model.Nnc*model.Ndp; %number of data per sample (nuclides*depths)
model.age = 3.0; %Max time (Myr)
model.z0 = 20; %Max depth
model.Temp = 1.0;
model.Mmp = 2; %number of generic model parameters
model.data{1}.density=2.65;

%load and prepare d18O data
load('d18Ocurves.mat','Age','d18O_4ky');
I = find(Age < model.age*1e6);
model.dO = d18O_4ky(I);
model.dOt = Age(I);

%% Loop scenarios to find exhumation histories
for i=1:Nscn
    T1(i)=Tdgla(i);
    T2(i) = T1(i) + dT2(i);
    Z2(i) = Z1(i) + 10^dZ2(i);
    
    T3(i) = T2(i) + dT3(i);
    Z3(i) = Z2(i) + 10^dZ3(i);
    
    T4(i) = model.age; %this requires age > Tdg+dT2+dT3
    Z4(i) = Z3(i) + (model.age - T3(i))*10^E0(i);
end

%% Prep figure layout
tiledlayout(2,4,"TileSpacing","compact","Padding","compact")

nexttile(1), hold on, box on %Exhumation histories
for i=1:Nscn
    up = [d18Op(i),Tdgla(i),Z1(i),dT2(i),dZ2(i),dT3(i),dZ3(i),E0(i)];
    syn_exhu_(up,lcols(i,:),mcols(i,:),1,lstyls{i},mstyls{i},i)
end
xlabel('Time before present (Ma)'), ylabel('Depth (m)')
set(gca,'ydir','reverse','FontSize',14)
ylim([0 3]), xlim([0 3])

nexttile(2,[1,3]), hold on, box on %Exposure histories
load d18Ocurves.mat
plot(Age/1e6,d18O_4ky,'k')
for i=1:Nscn
    line([0 max(xlim)],[d18Op(i) d18Op(i)],'Color',lcols(i,:),...
        'LineStyle',lstyls{i},'Linewidth',1)
    d=d18O_4ky; d(d18O_4ky<d18Op(i))=1; d(d18O_4ky>=d18Op(i))=2;
    dd=diff(d); Idd=find(abs(dd)>0);

    for j=1:floor(length(Idd)/2)
        patch([Age(Idd(2*j-1))/1e6 Age(Idd(2*j))/1e6 Age(Idd(2*j))/1e6...
            Age(Idd(2*j-1))/1e6],[5.2 5.2 5.28 5.28]+(i-1)*0.13,...
            brighten(mcols(i,:),0.9),'EdgeColor',lcols(i,:))
    end
    text(2.8,5.23+(i-1)*0.13,[num2str(round(sum(d(1:501)-1)/500*100,0)) '%'],...
    'color',lcols(i,:)) %Add burial since 500 ka in %
end
xlim([0 3]), ylim([2.75 5.5])
patch([min(xlim) max(xlim) max(xlim) min(xlim)],[5.08 5.08 5.15 5.15],'w','Edgecolor','k')
line([min(xlim) min(xlim)],[5.08 5.15],'color','w')
line([max(xlim) max(xlim)],[5.08 5.15],'color','w','linewidth',2)
xlabel('Time before present (Ma)'), ylabel('$\delta_{18}$O ($10^{-3}$)')
set(gca,'FontSize',14,'XTick',0:0.5:3,'Ytick',3:5,'ydir','reverse')

nexttile(5,[1,2]), hold on, box on %Two-nuclide diagram
banana_fig1(sample{1}.production)
text(5.75e5,5.6,'0.5 Myr')
text(5.75e5,4.43,'1.0 Myr')
text(5.75e5,3.5,'1.5 Myr')
xlabel('$^{10}$Be (at g$^{-1}$)'), ylabel('$^{26}$Al/$^{10}$Be')
set(gca,'FontSize',14,'xscale','log')
axis([5e5,5e7,3,8]);

nexttile(7), hold on, box on %10Be depth profile
xlabel('$^{10}$Be (at g$^{-1}$)'), ylabel('Depth (m)')
set(gca,'ydir','reverse','FontSize',14)

nexttile(8), hold on, box on %26Al/10Be ratio
xlabel('$^{26}$Al/$^{10}$Be'), ylabel('Depth (m)')
set(gca,'ydir','reverse','FontSize',14)


%% PLot pathways in two-nuclide diagram
nexttile(5)
dt=1e3; nt=ceil(model.age*10^6/dt); %number of timesteps
time = model.age*10^6*linspace(0,1,nt);
dOn = interp1(model.dOt,model.dO,time,'linear',model.dO(end));
for i=1:Nscn
    maxdepth=Z4(i)+1; nd=1e5; zd=maxdepth*linspace(0,1,nd); %depth vector
    NBe=zeros(nd,nt); NAl=zeros(nd,nt); %initiate Be10 and Al26 matrices
    erate=10^E0(i)/10^6; %change from m/Myr to m/yr
    [NBe(:,nt),NAl(:,nt)] = cosmoStepDepth(0,0,dt,0,...
        sample{1}.production,erate,zd,'ss'); %steady-state surface conc.
    mT = [0,T1(i),T2(i),T3(i),T4(i)]*1e6; %Times, Myr to yr
    mz = [0,Z1(i),Z2(i),Z3(i),Z4(i)]; %Depths, m
    burial = interp1(mT,mz,time); %Surface depths (m) at times in 'time'
    pfac = ones(nt,1); %controls exposure, modified below
    pfac(dOn > d18Op(i)) = 0; %No exposure when d18O values above threshold (d18Op)
    pfac(time < 25e3) = 0; %correct exposure around deglaciation: set no
                           %exposure last 25 kyr
    pfac(time < T1(i)*1e6) = 1; %then add exposure from time of deglaciation
    for kk=(nt-1):-1:1 %loop through time to get surface concentrations
        dt = time(kk+1)-time(kk); %timestep length
        dz = burial(kk+1)-burial(kk); %depth change
        erate = dz/dt; %erosion rate in time step
        pf = .5*(pfac(kk+1)+pfac(kk)); %exposure (production factor)
        [~,idx]=min(abs(zd-dz)); %Note: approximation for illustration only
        NBeinh=[NBe(idx:end,kk+1); NBe(1:idx-1,kk+1)]; %Inherited Be
        NAlinh=[NAl(idx:end,kk+1); NAl(1:idx-1,kk+1)]; %Inherited Al
        [NBe(:,kk),NAl(:,kk)] = cosmoStepDepth(NBeinh,NAlinh,...
            dt,pf,sample{1}.production,erate,zd,'eb'); %New surface conc.
    end
    plot(NBe(1,:),NAl(1,:)./NBe(1,:),'Color',lcols(i,:),'LineWidth',1.5,...
        'LineStyle',lstyls{i}) %plot surface concentration over time
end

%% Calculate and plot banana end points and depth profiles 
for i=1:Nscn
    up = [d18Op(i),Tdgla(i),Z1(i),dT2(i),dZ2(i),dT3(i),dZ3(i),E0(i)];
    [gm] = forward_bedrockvJ4(up,model,CNprop);
    nexttile(5)
    plot(gm(1),gm(model.Ndp+1)/gm(1),mstyls{i},"MarkerFaceColor",brighten(mcols(i,:),.9),...
        'MarkerEdgeColor',lcols(i,:),'MarkerSize',10,'LineWidth',1)
    nexttile(7)
    plot(gm(1:model.Ndp),model.data{1}.depths,mstyls{i},"MarkerFaceColor",...
        brighten(mcols(i,:),.9),'MarkerEdgeColor',lcols(i,:),...
        "LineStyle",lstyls{i})
    nexttile(8)
    plot(gm(model.Ndp+1:end)./gm(1:model.Ndp),model.data{1}.depths,...
        mstyls{i},"MarkerFaceColor",brighten(mcols(i,:),.9),...
        'MarkerEdgeColor',lcols(i,:),"LineStyle",lstyls{i})
end
nexttile(1), legend('Scenario A','Scenario B','Location','southwest')
nexttile(7), set(gca,'xscale','log')

%% Finish up
nexttile(1)
text(max(xlim)-0.1,min(ylim),'a','FontSize',20,'Horiz','right','Vert','top')
nexttile(2)
text(min(xlim)+0.02,min(ylim)+0.04,'b','FontSize',20,'Horiz','left','Vert','top')
nexttile(5)
text(max(xlim)-3e6,max(ylim),'c','FontSize',20,'Horiz','right','Vert','top')
nexttile(7)
text(min(xlim)+3e4,min(ylim),'d','FontSize',20,'Horiz','left','Vert','top')
nexttile(8)
text(max(xlim)-0.1,min(ylim),'e','FontSize',20,'Horiz','right','Vert','top')

set(gcf,'units','normalized','position',[.1,.3,.75,.75]);

%% Set interpreter to default
set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');