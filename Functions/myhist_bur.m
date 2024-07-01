function myhist_bur(x0,dx,ymin,ymax,ny,vals,col,smval,alpha,ls,lw,bur)

yi = linspace(ymin,ymax,ny); dy=(ymax-ymin)/ny;
[N,~]=histcounts(vals,yi); %[N2,~]=histcounts(uval2,yi);[N3,~]=histcounts(uval3,yi);
[Ncorr,~]=histcounts(bur,yi);
N=smooth(N./Ncorr,smval);
Nmax=max(N);
N=N/Nmax; %N2=N2/Nmax; N3=N3/Nmax;
BinCenters=yi(2:end)-0.5*dy;
%patch([smooth(N2(N2>0));zeros(length(N2(N2>0)),1)],[BinCenters(N2>0) fliplr(BinCenters(N2>0))], .7*[1 .1 .1]) %smooth histogram patch
%patch([smooth(N3(N3>0));zeros(length(N3(N3>0)),1)],[BinCenters(N3>0) fliplr(BinCenters(N3>0))], .7*[.3 1 .3], 'Linewidth', 2) %smooth histogram patch

patch([x0-dx*N(N>0);x0+flipud(dx*N(N>0))],[BinCenters(N>0) ...
    fliplr(BinCenters(N>0))], col,'EdgeColor','none','FaceAlpha',alpha) %smooth histogram patch
%    patch([x0-f(:);x0*ones(length(f),1)],[yi(:);flipud(yi(:))],col,'linestyle','none');

if alpha == 0
    line([x0-dx*N(N>0);x0+flipud(dx*N(N>0))],[BinCenters(N>0) ...
    fliplr(BinCenters(N>0))],'LineWidth',lw,'Color',col,'LineStyle',ls)
end