function syn_exhu_(up,cg,cig,lw,ls,ms,scn)
% Modified looks, JLA, 19.02.21
%colors
% cg = [0,.3,.2];
% cig = [.1,.4,.3];
xlim = get(gca,'xlim');

%depth points
zp=[0.0,up(3),up(3)+10^(up(5)),up(3)+10^(up(5))+10^(up(7)),...
    up(3)+10^(up(5))+10^(up(7))+(3 - (up(2)+up(4)+up(6)))*10^(up(8))];
tp=[0.0,up(2),(up(2)+up(4)),(up(2)+up(4)+up(6)),3];
% zend = zp(end) + (xlim(2)-tp(end))*(zp(end)-zp(end-1))/(tp(end)-tp(end-1));

%draw line segment
% line(tp,zp,'color',cg,'linewidth',lw,'linestyle',ls,'handlevisibility','off');
%         
% %draw points
mz = 10;
% plot(tp(2:(end-1)),zp(2:(end-1)),ms,'markeredgecolor',cg,'markersize',mz,...
%     'markerfacecolor',brighten(cig,.9),'LineWidth',1,'DisplayName',...
%     ['Scenario ' num2str(scn)]);

plot(tp,zp,'o-',"Color",cg,'markeredgecolor',cg,'markersize',mz,...
    'markerfacecolor',brighten(cig,.9),'LineWidth',lw,'DisplayName',...
    ['Scenario ' num2str(scn)],'linestyle',ls,'Marker',ms)

% pos = get(gca,'position');
% xlim = get(gca,'xlim');
% ylim = get(gca,'ylim');


% for j=3:length(tp)-1
%     xtp(j-1) = pos(1) + pos(3)*(tp(j)-xlim(1))/(xlim(2)-xlim(1));
%     ytp(j-1) = pos(2) + pos(4) + 0.02;
%     tlabel{j-1} = ['$T_',num2str(j-1),'$'];
%     xzp(j-1) = pos(1) + pos(3) + 0.02;
%     yzp(j-1) = pos(2) + pos(4)*(ylim(2)-zp(j))/(ylim(2)-ylim(1));
%     zlabel{j-1} = ['$Z_',num2str(j-1),'$'];
% end

