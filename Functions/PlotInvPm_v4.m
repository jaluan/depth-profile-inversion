function PlotInvPm_v4(fnr,tile,nScn)

%% Plot inversion of surface sample
load(['models/Syn_v4_scn',num2str(fnr),'_truesurf_n.mat']); %surface sample only

nexttile(tile) %Ice burial
load('burialTime_500ka.mat') %conversion of d18O to burial time in % within last 500 ka
dx=0.9; Nkernel = 500;
bvalt=interp1(dO_i,bt,model.synpmval(1)); %true burial value (% since 1 Ma)
line([0 1],[bvalt bvalt],'Color',[.8,.3,.3],'LineWidth',2)% true burial value
hold on
[~,~,~,bval,~,~]=IceHist(model,dO_i,bt);

% Correcting for non-linear relation between dO_i and bt
a = 3.5; b = 5.0; %range
rdO18 = (b-a).*rand(1e6,1) + a; %random d18O values within range
bur=interp1(dO_i,bt,rdO18); %random burial values within range
myhist_bur(0,dx,0,100,Nkernel,bval,'k',5,.2,':',2,bur);%histogram of modelled results (corrected)
set(gca,'FontSize',16); box on
xlim([0 1]), ylim([0 100])
title(['Scenario ', num2str(fnr)])

nexttile(tile+nScn) %Time of 1 m erosion
TE1=0.5; %Ma
ZT1=1; %m
[E1t,TZ1t] = trueModel(model.synpmval,model.age,TE1,ZT1);  % true value
line([0 1],[TZ1t TZ1t],'Color',[.8,.3,.3],'LineWidth',2)
hold on
[E1,TZ1] = interpModel(model,TE1,ZT1); %model values
myhist_(0,dx,0,3,Nkernel,TZ1,'k',5,.2,':',2);%histogram of modelled results
set(gca,'FontSize',16); box on
xlim([0 1]), ylim([0 3])

nexttile(tile+2*nScn) %Erosion since 500 ka
Emax=5;
line([0 1],[E1t E1t],'Color',[.8,.3,.3],'LineWidth',2) % true value
hold on
myhist_(0,dx,0,Emax,Nkernel,E1,'k',5,.2,':',2);%histogram of modelled results
xlim([0 1]), ylim([0 Emax])
set(gca,'FontSize',16); box on

%% True concentration inversions depth profiles
load(['models/Syn_v4_scn',num2str(fnr),'_true_n.mat']); % depth profile inversion, true values
nexttile(tile) %Ice burial
[~,~,~,bval,~,~]=IceHist(model,dO_i,bt);
myhist_bur(0,dx,0,100,Nkernel,bval,[.3,.3,.6],5,0,':',3,bur);%histogram of modelled results

nexttile(tile+nScn) %Time of 1 m erosion
[E1,TZ1] = interpModel(model,TE1,ZT1); %model values
myhist_(0,dx,0,3,Nkernel,TZ1,[.3,.3,.6],5,0,':',3);%histogram of modelled results

nexttile(tile+2*nScn) %Erosion since 500 ka
myhist_(0,dx,0,Emax,Nkernel,E1,[.3,.3,.6],5,0,':',3);%histogram of modelled results