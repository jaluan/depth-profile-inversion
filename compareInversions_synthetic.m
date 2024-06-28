% This script creates a figure to compare inversions based on surface
% sample vs. depth profile for synthetic and real data

close all;

set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

addpath('Functions')

% CNprop = getCNprop; % Cosmogenic halflives etc.
% load colormaps.mat

%% load data

scn=1; %Change scenario here

samplename=['Scenario ', num2str(scn)];
load(['models/Syn_v4_scn',num2str(scn),'_truesurf_n.mat']); %surface sample only
model1=model;
load(['models/Syn_v4_scn',num2str(scn),'_true_n.mat']); %depth points incl. surface sample, true values
model2=model;
clear model

%% Exhumation plots
Maxz = 5; %y-axis max depth
Maxt = model1.age;
Nz = 100;
Nt = 200;
zint = linspace(0,Maxz,Nz);
tint = linspace(0,Maxt,Nt);
[tbin,zbin]=meshgrid(tint,zint);
cvals = linspace(-4,0,100);

% Plot exhumation paths for surface sample
ax1=axes(); set(gca,'ydir','reverse'); 
% set(gca,'Color',[.7 .7 .7]);%set(gca,'visible','off');
hold on; box on; set(gca,'layer','top');
map = .8*[1 1 1;1 1 1;1 1 1]; %gray; 
map(1,:) = [0,0,0]; %[1,1,1],  %colormap
colormap(ax1,map);
histgrid=hgrid(tbin,zint,tint,model1);
% [C1,h1] = contourf(ax1,tbin,zbin,log10(histgrid+1e-6),cvals,'linestyle','none');
[C1,h1] = contourf(ax1,tbin,zbin,log10(histgrid+1e-6),[cvals(1),cvals(1)],'linestyle','none');
% [C2,h2] = contour(ax1,tbin,zbin,log10(histgrid+1e-6),[cvals(1),cvals(1)],'color','k');

% Plot exhumation paths for depth profile
ax2=axes(); set(gca,'ydir','reverse'); 
set(ax2,'Color','none');
hold on; box on; set(gca,'layer','top');
% create a color map ranging from blue to light pink
    length = 255;
    blue = [.3 .3 .6];
    pink = [.8 .55 .7]; %[1 .75 .8];
    map = [linspace(blue(1),pink(1),length)', ...
        linspace(blue(2),pink(2),length)', linspace(blue(3),pink(3),length)'];
    colormap(ax2,map);
% map = bone; map(1,:) = [0,0,0]; %[1,1,1]; % bone colormap
% map = [.3 .3 .6;.3 .3 .6;.3 .3 .6]; map(1,:) = [0,0,0]; %[1,1,1]; %colormap
colormap(ax2,map);
histgrid=hgrid(tbin,zint,tint,model2);
contourf(ax2,tbin,zbin,log10(histgrid+1e-6),cvals,'linestyle','none'); %colormap
% contourf(ax2,tbin,zbin,log10(histgrid+1e-6),[cvals(1),cvals(1)],'linestyle','none');

% % Plot exhumation paths for depth profile with extra depth points
% ax3=axes(); set(gca,'ydir','reverse'); 
% set(ax3,'Color','none');
% hold on; box on; set(gca,'layer','top');
% map = pink; map(1,:) = [0,0,0];%[1,1,1] %colormap
% colormap(ax3,map);
% histgrid=hgrid(tbin,zint,tint,model3);
% contourf(ax3,tbin,zbin,log10(histgrid+1e-6),cvals,'linestyle','none');

% Plot true exhumation path
ax4=axes(); set(gca,'ydir','reverse'); 
set(ax4,'Color','none');
hold on; box on; set(gca,'layer','top');
syn_exhu(model1.synpmval)
ylim([0 5])

% Draw stipled line at summary figure cuts
plot([.5 .5],[0 5],'k','LineStyle','--','Linewidth',2)
% text(.8,1.8,'Erosion since 500 ka ','Rotation',90,'Interpreter','latex','Fontsize',20);
text(0.57,4.7,'Erosion since 500 ka ','Rotation',90,'Interpreter','latex','Fontsize',30);
plot([0 3],[1 1],'k','LineStyle','--','Linewidth',2)
text(2.1,0.85,'Time of 1 m erosion','Interpreter','latex','Fontsize',30);

set([ax1,ax2,ax4],'Position',[0.1,0.1,0.85,0.85]);
set(ax1,'xtick',[]); set(ax1,'ytick',[]);
set(ax2,'xtick',[]); set(ax2,'ytick',[]);
% set(ax3,'xtick',[]); set(ax3,'ytick',[]);
set(ax4,'xtick',0:1:5,'fontsize',35); set(ax4,'ytick',0:1:5,'fontsize',30);

xlabel('Time (Ma)','fontsize',30);
ylabel('Depth below surface (m)','fontsize',30);
% title(samplename,'FontSize',30)
text(min(xlim)+(max(xlim)-min(xlim))/2,min(ylim),samplename,'FontSize',30,'Horiz','center','Vert','top')

%% Draw ice histories as violin plots in inset
dx=0.3; Nkernel = 500;

% ax5 = axes('position',[0.79,0.12,0.15,0.25]); %d18O
% hold on; box on;
% set(gca,'Color',1*[1 1 1]) %'none'
% set(gca,'xtick',[]); set(gca,'ytick',[3.5;5]);
% set(gca,'yticklabel',[3.5;5],'fontsize',20,'YColor',.1*[1 1 1]);
% axis([0,3,3.5,5])

ax6 = axes('position',[0.67,0.15,0.25,0.39]); %ice burial
hold on; box on;
set(gca,'Color',1*[1 1 1]) %'none'
set(gca,'xtick',[]);
set(gca,'ytick',[0;100]);
set(gca,'yticklabel',[0;100],'fontsize',24); % set(gca,'ydir','reverse')
axis([0,3,0,100]);

load('burialTime_500ka.mat') %conversion of d18O to burial time in % within last 500 ka

% Plot true ice history
dO=model1.synpmval(1); bO=interp1(dO_i,bt,dO);
% axes(ax5), line([0 3],[dO dO],'Color',[.8,.3,.3],'LineWidth',4)
axes(ax6), line([0 3],[bO bO],'Color',[.8,.3,.3],'LineWidth',4)

% Plot ice history for surface sample
[~,~,~,bval,~,~]=IceHist(model1,dO_i,bt);
% axes(ax5), myhist(1,dx,3,5,Nkernel,uval,.8*[1 1 1],1);
axes(ax6), myhist(1,dx,0,100,Nkernel/5,bval,.8*[1 1 1],3);
% plot(1,max(bval),'o','Color',.6*[1 1 1],'MarkerSize',10,'LineWidth',2)
% plot(1,min(bval),'o','Color',.6*[1 1 1],'MarkerSize',10,'LineWidth',2)

% Plot ice history for depth profile with 5 data points
[uval,~,~,bval,~,~]=IceHist(model2,dO_i,bt);
% axes(ax5), myhist(2,dx,3,5,Nkernel,uval,[.3 .3 .6],1);
axes(ax6), myhist(2,dx,0,100,Nkernel/5,bval,[.3 .3 .6],3);
% plot(2,max(bval),'o','Color',[.3 .3 .6],'MarkerSize',10,'LineWidth',2)
% plot(2,min(bval),'o','Color',[.3 .3 .6],'MarkerSize',10,'LineWidth',2)

% % Plot ice history for depth profile with 10 data points
% [uval,uval2,uval3,bval,bval2,bval3]=IceHist(model3,dO_i,bt);
% axes(ax5), myhist(3,dx,3,5,Nkernel,uval,col3b,1);
% axes(ax6), myhist(3,dx,0,100,Nkernel/5,bval,col3b,5);

% Add labels
% axes(ax5),
% text(-.33,4,'$\delta^{18}$O$_{th}$','Interpreter','latex','Rotation',90,'Fontsize',30);
axes(ax6),
text(-.33,15,'Ice burial ($\%$)','Interpreter','latex','Rotation',90,'Fontsize',32);

%% Draw legend
% text(1.7,3.72,'Surface sample, [$^{10}$Be, $^{26}$Al]','Interpreter','latex','Fontsize',20);
% text(1.5,4.1,'$\bf{2 \ m \ depth \ profiles, \ [^{10}Be, ^{26}Al]:}$','Interpreter','latex','Fontsize',20);
% text(1.7,4.37,'+4 depth points','Interpreter','latex','Fontsize',20);
% text(1.7,4.72,'+9 depth points','Interpreter','latex','Fontsize',20);
% patch([1.5 1.5 1.6 1.6],[3.6 3.8 3.8 3.6],[.3 .3 .41]) %blue patch
% patch([1.5 1.5 1.6 1.6],[3.6 3.8 3.8 3.6],[.3 .3 .41],'FaceColor','none','EdgeColor','k','LineStyle','--','Linewidth',2) %stipled black line
% patch([1.5 1.5 1.6 1.6],[4.25 4.45 4.45 4.25],[.3 .65 .4]) % green patch
% line([1.49 1.61],[4.35 4.35],'Color',[0,.3,.2],'LineWidth',2) %green line
% plot(1.55,4.35,'o','MarkerFaceColor',[.1,.4,.3],'MarkerEdgeColor','k','MarkerSize',8) %green marker
% patch([1.5 1.5 1.6 1.6],[4.6 4.8 4.8 4.6],[.67 .44 .44]) %red patch
% line([1.49 1.61],[4.7 4.7],'Color',[.4,.1,.1],'LineWidth',2) %red line
% plot(1.55,4.7,'o','MarkerFaceColor',[.8,.4,.3],'MarkerEdgeColor','k','MarkerSize',8) %red marker

set(gcf,'units','normalized','position',[.1,.3,.8,.8]);
% mname = ['models/SyntheticInversion'];
% print(mname,'-dpdf','-bestfit');
% print(mname,'-dpng','-r300');

set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');