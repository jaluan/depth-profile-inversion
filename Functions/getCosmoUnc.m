function [Nunc] = getCosmoUnc(N,W,nuclide,pflag)

% Input
% N: Concentration of nuclide in at/g
% W: Weight of dissolved sample (quartz) in g
% nuclide: 'Be' or 'Al'
% pflag: 1=plot, else no plot

% Output
% Nunc: anticipated internal AMS measurement uncertainty in %

% JLA, Feb 13th 2020

Nav = 6.02214086*10^23; %Avogadros number, mol-1

[rows,cols]=size(N);
if cols>1
    N=reshape(N,rows*cols,1);
end

switch nuclide
    case 'Be10'
        N10=N*W; % Be10 atoms
        Mspike=250*10^-6; % Be9 spike, [g]
        MM9=9.012182; % Molar mass Be9, [g mol-1]
        N9=Mspike/MM9*Nav; % Be9 atoms
        Beratio=N10/N9; 
        
        [num,~,~]=xlsread('data/N10all.xlsx');
        I=num(:,6)>=20 & num(:,4)>=1; %only look at samples where runtime > 20 min, and total charge > 1 mC
        error=num(I,2);
        ratio=num(I,1)*10^-12;
        [B,I] = sort(ratio);
        C=error(I); 
        %M = movmean(C,5); %old moving mean method, too jumpy
        M=fit(B,C,'power1'); %power fit

        Nunc=zeros(size(N10));
        
        for i=1:length(N10)
            if Beratio(i)>=max(ratio)
             Nunc(i)=1.1;
            elseif Beratio(i)<=min(ratio)
                   Nunc(i)=200;
            else
                Nunc(i)=interp1(B,M(B),Beratio(i));
            end
        
            disp(['N10 uncertainty: ' num2str(round(Nunc(i),2)) '%']);
        
            if pflag == 1
              col=[0.494 0.184 0.556];
              scatter(ratio,error,20,'filled','MarkerFaceColor',col), hold on
              plot(B,M(B),'color',brighten(col,-0.6),'LineWidth',2,'LineStyle','-.')
              plot(Beratio(i),Nunc(i),'ro','MarkerSize',20,'HandleVisibility','off')
            end
        end
        
    case 'Al26'
        N26=N*W; % Al26 atoms
        
        Mspike=1200*10^-6; % Al27 spike, [g]
        MM27=728.501542; % Molar mass Al27, [g mol-1]
        N27=Mspike/MM27*Nav; % Al27 atoms
        Alratio=N26/N27; 
        
        [num,~,~]=xlsread('data/N26all.xlsx');
        I=num(:,6)>=30; %only look at samples where runtime > 30 min
        error=num(I,2);
        ratio=num(I,1)*10^-12;
        [B,I] = sort(ratio);
        C=error(I); 
        %M = movmean(C,5); %old moving mean method, too jumpy
        M=fit(B(B>1E-14),C(B>1E-14),'power1'); %power fit to conc>1e-14

        Nunc=zeros(size(N26));
        
        for i=1:length(N26)
            if Alratio(i)>=max(ratio)
             Nunc(i)=0.7;
            elseif Alratio(i)<=min(ratio)
                   Nunc(i)=200;
            else
                Nunc(i)=interp1(B,M(B),Alratio(i));
            end
        
            disp(['N26 uncertainty: ' num2str(round(Nunc(i),2)) '%']);
        
            if pflag == 1
              col=[0.929 0.694 0.125];
              scatter(ratio,error,20,'filled','MarkerFaceColor',col,'Marker','diamond'), hold on
              plot(B,M(B),'color',brighten(col,-0.6),'LineWidth',2)
              plot(Alratio(i),Nunc(i),'ko','MarkerSize',20,'HandleVisibility','off')
            end
        end
        
    otherwise
        Nunc=100;
        disp('Warning: This nuclide is not implemented yet')
end

if cols>1
    Nunc=reshape(Nunc,rows,cols);
end
