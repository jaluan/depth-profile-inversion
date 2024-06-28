function makereport_syn_dp(snr,fnr)

close all;
set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

addpath('Functions/export_fig')
% Cosmogenic properties
CNprop = getCNprop;


%load MC results
str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

mname = ['models/Syn_v4_scn',num2str(fnr),'_true_n.mat'];  dpname = ['Synv4_scn',num2str(fnr)];
load(mname);

col1 = 0.9*[1,1,1];
col2 = 0.6*[1,1,1];
col3 = 0.3*[1,1,1];

map = colormap;
nc = length(map);
xc = linspace(0,1,nc);
rc = map(:,1);
gc = map(:,2);
bc = map(:,3);
for i=1:model.Nwalk
    ri = interp1(xc,rc,(i-1)/model.Nwalk);
    gi = interp1(xc,gc,(i-1)/model.Nwalk);
    bi = interp1(xc,bc,(i-1)/model.Nwalk);
    wcol(i,:) = [ri,gi,bi];
end

%plot margins
lpm = 0.1;
ddm = 0.02;

%******* figures showing walkers *************

set(gcf,'units','normalized','position',[.1,.3,.4,.6]);
set(gca,'position',[0,0,1,1],'visible','off');
set(gca,'xlim',[0,1],'ylim',[0,1]);
set(gcf,'Name','Walker information');
ax0 = gca;

[np,nn] = numSubplots(model.Nmp);

%plot models params
for i = 1:model.Nmp
    subplot(np(1),np(2),i); hold on; box on; grid on;
    xlabel('model nr.');
    ylabel(model.mp{i}.name);
    set(gca,'xlim',[0,length(model.walker{1}.status)],'ylim',[model.mp{i}.vmin,model.mp{i}.vmax]); %
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == -2); %rejected during burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col1);
        I = find(model.walker{nw}.status == 0); %rejected after burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col2);
        I = find(model.walker{nw}.status == -1); %accepted during burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col3);
        I = find(model.walker{nw}.status == 1); %accepted after burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',wcol(nw,:));
    end
    plot([0,length(model.walker{1}.status)],[model.synpmval(i),model.synpmval(i)],'r-','LineWidth',2); %plot synthetic pm values
end

print('temp3.pdf','-dpdf','-fillpage');

%********************** Parameter distributions ************************
figure()
set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
set(gcf,'Name','Parameter distributions');

[np,nn] = numSubplots(model.Nmp);

%loop parameters
for i=1:model.Nmp
    
    subplot(np(1),np(2),i); 
    hold on; box on; grid on;
    ylabel('Frequency');
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
                    model.mp{i}.vmax]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status > -1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval);
    end
    line(xi,f,'color','k','linewidth',2);
    
    ylims = get(gca,'ylim');
    plot([model.synpmval(i),model.synpmval(i)],[ylims(1),ylims(2)],'r-','LineWidth',2); %plot synthetic pm values

end

print('temp2.pdf','-dpdf','-fillpage');

%21Ne figures
if model.data{1,1}.Nnuc==3 %if 21Ne in model
    % initiate figure
    figure()
    set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
    set(gcf,'Name','21Ne figures');
    ax0 = gca;
    set(ax0,'position',[0,0,1,1]);
    set(gca,'visible','off');

    dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;
    
    for ns = 1:model.Nsnr %loop samples

        axes('position',[lpm+(ns-1)*(ddm+dpx),0.7,dpx,0.25]);
        hold on; box on; grid on;
        title(model.data{ns}.name);

        n0 = (ns-1)*model.Nds; 
        Ndp = model.Ndp;

        % retrieve modelled nuclide distributions
        Be = [];
        Ne = [];
        Al = [];
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status == 1);
            Be = [Be(:,:);model.walker{nw}.gm(n0+1:Ndp+n0,I)'];
            Al = [Al(:,:);model.walker{nw}.gm(n0+Ndp+1:n0+2*Ndp,I)'];
            Ne = [Ne(:,:);model.walker{nw}.gm(n0+2*Ndp+1:n0+3*Ndp,I)'];
        end
        BeNe = Be./Ne; %calculate ratios
        AlNe = Al./Ne;
        Neint = linspace(0.5*min(min(Ne)),max(max(Ne)),40)';
        BeNeint = linspace(0.5*min(min(BeNe)),1.5*max(max(BeNe)),40)';
        AlNeint = linspace(0.5*min(min(AlNe)),1.5*max(max(AlNe)),40)';

        %Ne-Be banana
        xlabel(['Normalized NNe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N10Be/N21Ne');
        else, set(gca,'yticklabel',[]);
        end
        
        
        for j=1:Ndp
            N = hist3([BeNe(:,j),Ne(:,j)],{BeNeint Neint});
            N = N/sum(N(:));

            %Normalize to SLHL using surface production rate (spallation+muons)
            nfac = CNprop.PNe0/(model.data{ns}.production.P21spal + ...
                model.data{ns}.production.P21_m1 + model.data{ns}.production.P21_m2);
            [X,Y] = meshgrid(Neint*nfac,BeNeint);
            contour(X,Y,N,40);
        end

        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r1021,model.data{ns}.dN21*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r1021,model.data{ns}.dr1021,'vertical','.k');
        
        %Ne-Al banana
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.075,dpx,0.25]);
        hold on; box on; grid on;
        
        xlabel(['Normalized NNe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N26Al/N21Ne');
        else, set(gca,'yticklabel',[]);
        end
        
        for j=1:Ndp %loop depthpoints
            N = hist3([AlNe(:,j),Ne(:,j)],{AlNeint Neint});
            N = N/sum(N(:));

            %Normalize to SLHL using surface production rate (spallation+muons)
            nfac = CNprop.PNe0/(model.data{ns}.production.P21spal + ...
                model.data{ns}.production.P21_m1 + model.data{ns}.production.P21_m2);

            [X,Y] = meshgrid(Neint*nfac,AlNeint);
            contour(X,Y,N,40);
        end
        
        model.data{ns}.r2621 = model.data{ns}.N26./model.data{ns}.N21;
        model.data{ns}.dr2621 = model.data{ns}.r2621.*sqrt((model.data{ns}.dN26./...
        model.data{ns}.N26).^2+(model.data{ns}.dN21./model.data{ns}.N21).^2);
    
        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r2621,model.data{ns}.dN21*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r2621,model.data{ns}.dr2621,'vertical','.k');
        
        % Histogram/kernel density with modelled concentrations
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.4,dpx,0.2]);
        hold on; box on; grid on;
        xlabel('N21 (atoms/g)'); ylabel ('Probability')
        histogram(Ne,'Normalization','probability')

        % Measured data on top
        errorbar(model.data{ns}.N21,.015*ones(Ndp,1),model.data{ns}.dN21,...
        'horizontal','.k','Linewidth',1.5,'Color',[.7 .2 .2]);
    
    end
    print('temp0.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below
end

%******* report figure *************
figure;
set(gcf,'papertype','a4');
set(gcf,'units','centimeters','position',[5,5,21,29.7]);
set(gcf,'Name','Report');

ax0 = gca;
set(ax0,'position',[0,0,1,1]);
set(gca,'visible','off');
text(0.05,0.97,['File: ',mname],'HorizontalAlignment','left','fontsize',14);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

for ns = 1:model.Nsnr

    axes('position',[lpm+(ns-1)*(ddm+dpx),0.675,dpx,0.25]);
    hold on; box on; grid on;
%     set(gca,'ylim',[3,7.5]);
    
    title(model.data{ns}.name);
    
    n0 = (ns-1)*model.Nds;
    Ndp = model.Ndp;
    
    Be = [];
    Al = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status > -1);
        Be = [Be(:,:);model.walker{nw}.gm(n0+1:Ndp+n0,I)'];
        Al = [Al(:,:);model.walker{nw}.gm(n0+Ndp+1:n0+2*Ndp,I)'];
    end
    
    AlBe = Al./Be;
    
    Beint = linspace(0.5*min(min(Be)),max(max(Be)),40)';
    Alint = linspace(min(min(Al)),max(max(Al)),40)';
    AlBeint = linspace(0.5*min(min(AlBe)),1.5*max(max(AlBe)),40)';
    
    %set(gca,'xlim',[0.5*min(Be(:)),1.5*max(Be(:))]);
    
    %banana(CNprop);
    xlabel(['NBe (atoms/g), sample ',num2str(ns)]);
    if (ns == 1), ylabel(['NAl/NBe']);
    else, set(gca,'yticklabel',[]);
    end
    
    
    for j=1:Ndp
        N = hist3([AlBe(:,j),Be(:,j)],{AlBeint Beint});
        N = N/sum(N(:));

        nfac = CNprop.PBe0/(model.data{ns}.production.P10spal + model.data{ns}.production.P10_m1+ model.data{ns}.production.P10_m2);
        [X,Y] = meshgrid(Beint*nfac,AlBeint);
        contour(X,Y,N,40);
    end
    
    r2610true=model.data{ns}.N26true./model.data{ns}.N10true;
    dr2610true=r2610true.*sqrt((model.data{ns}.dN26true./model.data{ns}.N26true).^2 ...
        +(model.data{ns}.dN10true./model.data{ns}.N10true).^2);
    errorbar(model.data{ns}.N10true*nfac,r2610true,model.data{ns}.dN10true*nfac,'horizontal','.k');
    errorbar(model.data{ns}.N10true*nfac,r2610true,dr2610true,'vertical','.k');
  
end

%****** glaciation history *********

dpx = (1-2*lpm-ddm)/2;

%loop parameters
for i=1:2
    
    axes('position',[lpm+(i-1)*(ddm+dpx),0.37,dpx,0.25]);
    hold on; box on; grid on;
    
    if (i == 1), ylabel('Frequency'); end
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
                    model.mp{i}.vmax]);
    set(gca,'yticklabel',[]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status > -1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval);
    end
    line(xi,f,'color','k','linewidth',2);
    ylims = get(gca,'ylim');
    plot([model.synpmval(i),model.synpmval(i)],[ylims(1),ylims(2)],'r-','LineWidth',2); %plot synthetic pm values
    
end


%exhumation plots
Maxz = model.z0;
Maxt = model.age;
z0 = model.z0;
Nz = 100;
Nt = 200;
zint = linspace(0,Maxz,Nz);
tint = linspace(0,Maxt,Nt);
[tbin,zbin]=meshgrid(tint,zint);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

map = colormap;
map(1,:) = [1,1,1];
colormap(map);

for ns = 1:model.Nsnr

    %subplot(np(1),np(2),ns); 
    axes('position',[lpm+(ns-1)*(ddm+dpx),0.05,dpx,0.25]);
    hold on; box on; set(gca,'layer','top');
    
    set(gca,'ydir','reverse');
    xlabel('Time (Ma)');
    if (ns == 1) 
        ylabel('Burial depth');
    else
        set(gca,'yticklabel',[]);
    end
    
    histgrid = zeros(size(tbin));
    dzi = zint(2)-zint(1);
    title(model.data{ns}.name);
    
    n0 = model.Mmp + (ns-1)*model.Nsmp;
    
    for nw=1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        for i=1:length(I)
            
           T1 = model.walker{nw}.u(2,I(i));
            z1 = model.walker{nw}.u(n0+1,I(i));
            dT2 = model.walker{nw}.u(n0+2,I(i));
            dz2 = 10^model.walker{nw}.u(n0+3,I(i));
            dT3 = model.walker{nw}.u(n0+4,I(i));
            dz3 = 10^model.walker{nw}.u(n0+5,I(i));
            E4 = 10^model.walker{nw}.u(n0+6,I(i));
            
            T2 = T1 + dT2;
            z2 = z1 + dz2;
            T3 = T2 + dT3;
            z3 = z2 + dz3;
            T4 = model.age; %this requires age > Tdg+dT2+dT3
            z4 = z3 + (model.age - T3)*E4;
        
            Tm = [0,T1,T2,T3,T4];
            zm = [0,z1,z2,z3,z4];
            
            zinterp = interp1(Tm,zm,tint,'linear',2*z0);        
            izs = 1+floor(zinterp/dzi); %bin index at each tsfine
            for itsfine=1:Nt
                if izs(itsfine)<=Nz %&& izs(itsfine)>0) %JLA added second clause 11.03.20
                    histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                end
            end
        end
    end
    
    
    [C,h] = contourf(tbin,zbin,(histgrid+0.5).^0.25,50, ...
                     'linestyle','none');
    
                
end
syn_exhu(model.synpmval)
ylim([0 5])

mname = ['models/reports/Report_sample_',str,'_',dpname,'.pdf'];

print(mname,'-dpdf','-fillpage');

if model.data{1,1}.Nnuc==2
    append_pdfs(mname, 'temp2.pdf', 'temp3.pdf')
    delete('temp2.pdf','temp3.pdf')
elseif model.data{1,1}.Nnuc==3
    append_pdfs(mname,'temp0.pdf', 'temp2.pdf', 'temp3.pdf')
    delete('temp0.pdf','temp2.pdf','temp3.pdf')
end

set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');