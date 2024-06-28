set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

close all;
addpath('Functions')

nScn=8; %number of scenarios
tiledlayout(3,nScn,"TileSpacing","compact",'Padding','compact')

% Create Fig. 4a
PlotInvPm_v4(1,1,nScn)
PlotInvPm_v4(2,2,nScn)
PlotInvPm_v4(3,3,nScn)
PlotInvPm_v4(4,4,nScn)
PlotInvPm_v4(9,5,nScn)
PlotInvPm_v4(10,6,nScn)
PlotInvPm_v4(11,7,nScn)
PlotInvPm_v4(12,8,nScn)

% % Create Fig. 4b
% PlotInvPm_v4(5,1,nScn)
% PlotInvPm_v4(6,2,nScn)
% PlotInvPm_v4(7,3,nScn)
% PlotInvPm_v4(8,4,nScn)
% PlotInvPm_v4(13,5,nScn)
% PlotInvPm_v4(14,6,nScn)
% PlotInvPm_v4(15,7,nScn)
% PlotInvPm_v4(16,8,nScn)

nexttile(1) %Ice burial
ylabel('Ice burial since 500 ka ($\%$)')

nexttile(nScn+1) %Time of 1 m erosion
hyl1 = ylabel('Time of 1 m erosion (Ma)');

nexttile(2*nScn+1) %Erosion since 500 ka
hyl2=ylabel('Erosion since 500 ka (m)');

for i=1:2*nScn
nexttile(i), set(gca,'Xtick',[])
end
for i=2:nScn
nexttile(i), set(gca,'Ytick',[])
end
for i=nScn+2:2*nScn
nexttile(i), set(gca,'Ytick',[])
end
for i=2*nScn+2:3*nScn
nexttile(i), set(gca,'Ytick',[])
end
%% Finish up
% set(gcf,'units','normalized','position',[.05,.1,nScn/15,.6]);
set(gcf,'units','normalized','position',[0.1,0.1,0.65,0.75]);
hyl1.Position(1)=hyl1.Position(1)-0.18; %adjust ylabel position
hyl2.Position(1)=hyl2.Position(1)-0.16; %adjust ylabel position

set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');