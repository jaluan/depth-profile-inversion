
function banana_fig1(prod)

% banana() produces a cosmo banana plot
% 
% DLE 11/10 2016
% JLA 10.03.23

% set(gca,'xscale','log');
% hold on; box on;
% axis([4e4,1e7,1,8]);


%************* Cosmo info **************
P10spal = prod.P10spal; %atoms/(g*yr)
P10m1 = prod.P10_m1;
P10m2 = prod.P10_m2;
P10m3 = prod.P10_m3;
P10tot = P10spal + P10m1 + P10m2 + P10m3;

P26spal = prod.P26spal; %atoms/(g*yr)
P26m1 = prod.P26_m1;
P26m2 = prod.P26_m2;
P26m3 = prod.P26_m3;
P26tot = P26spal + P26m1 + P26m2 + P26m3;

Lspal = prod.Lspal; % Attenuation lengths g/cm^2
P10Lm1 = prod.P10_Lm1;
P10Lm2 = prod.P10_Lm2;
P10Lm3 = prod.P10_Lm3;
P26Lm1 = prod.P26_Lm1;
P26Lm2 = prod.P26_Lm2;
P26Lm3 = prod.P26_Lm3;

TBe = 1.387e6;
TAl = 0.705e6;
lambda_Be = log(2)/TBe;
lambda_Al = log(2)/TAl;
rho = 2.65; %density g/cm3
            
%**************************************
time = linspace(0,1e7,1000)';
erates = [0,1e-7,1e-6,5e-6,1e-5,2e-5,1e-3];

NBe = zeros(length(time),length(erates));
NAl = zeros(length(time),length(erates));

for i=1:length(erates)

    %Full exposure with erosion
    fBe = lambda_Be + rho*erates(i)*100/Lspal;
    fAl = lambda_Al + rho*erates(i)*100/Lspal;
    NBe(:,i) = P10spal/fBe*(1-exp(-fBe*time));
    NAl(:,i) = P26spal/fAl*(1-exp(-fAl*time));

    fBe = lambda_Be + rho*erates(i)*100/P10Lm1;
    fAl = lambda_Al + rho*erates(i)*100/P26Lm1;
    NBe(:,i) = NBe(:,i) + P10m1/fBe*(1-exp(-fBe*time));
    NAl(:,i) = NAl(:,i) + P26m1/fAl*(1-exp(-fAl*time));

    fBe = lambda_Be + rho*erates(i)*100/P10Lm2;
    fAl = lambda_Al + rho*erates(i)*100/P26Lm2;
    NBe(:,i) = NBe(:,i) + P10m2/fBe*(1-exp(-fBe*time));
    NAl(:,i) = NAl(:,i) + P26m2/fAl*(1-exp(-fAl*time));

    fBe = lambda_Be + rho*erates(i)*100/P10Lm3;
    fAl = lambda_Al + rho*erates(i)*100/P26Lm3;
    NBe(:,i) = NBe(:,i) + P10m3/fBe*(1-exp(-fBe*time));
    NAl(:,i) = NAl(:,i) + P26m3/fAl*(1-exp(-fAl*time));
    
    
end
    
%No erosion    
line(NBe(:,1),NAl(:,1)./NBe(:,1),'color','k','linewidth',2);

% %with erosion
% for i=2:length(erates)
%     line(NBe(:,i),NAl(:,i)./NBe(:,i),'color','m','linestyle','--');
% end
% line(NBe(:,2),NAl(:,2)./NBe(:,2),'color','k','linestyle','--');
% line(NBe(:,3),NAl(:,3)./NBe(:,3),'color','k','linestyle','--');
% line(NBe(:,4),NAl(:,4)./NBe(:,4),'color','k','linestyle','--');
% line(NBe(:,5),NAl(:,5)./NBe(:,5),'color','k','linestyle','--');
% line(NBe(:,6),NAl(:,6)./NBe(:,6),'color','k','linestyle','--');


%infinite time with erosion
fero = logspace(-8,-2,200);
NBe_f = P10spal./(lambda_Be + rho*fero/Lspal)+...
    P10m1./(lambda_Be + rho*fero/P10Lm1)+...
    P10m2./(lambda_Be + rho*fero/P10Lm2)+...
    P10m3./(lambda_Be + rho*fero/P10Lm3);
NAl_f = P26spal./(lambda_Al + rho*fero/Lspal)+...
    P26m1./(lambda_Al + rho*fero/P26Lm1)+...
    P26m2./(lambda_Al + rho*fero/P26Lm2)+...
    P26m3./(lambda_Al + rho*fero/P26Lm3);
line(NBe_f,NAl_f./NBe_f,'color','k');

%burial isolines
tburial = [0.5e6,1e6,1.5e6,2e6,2.5e6,3e6];
NBe_b = zeros(length(time),length(tburial));
NAl_b = zeros(length(time),length(tburial));
for i=1:length(tburial)
   NBe_b(:,i) = NBe(:,1).*exp(-lambda_Be*tburial(i)); 
   NAl_b(:,i) = NAl(:,1).*exp(-lambda_Al*tburial(i)); 
   line(NBe_b(:,i),NAl_b(:,i)./NBe_b(:,i),'color','k','linestyle','--');
end

%burial paths
tb = linspace(0,3e6,200);
etime = [1e2,1e3,1e4,1e5,1e6,1e7]; %preburial exposure times
NBe_p = zeros(length(tb),length(etime));
NAl_p = zeros(length(tb),length(etime));
for i=1:length(etime)
    NBe_p(:,i) = P10tot/lambda_Be*(1-exp(-lambda_Be*etime(i)))*exp(-lambda_Be*tb);
    NAl_p(:,i) = P26tot/lambda_Al*(1-exp(-lambda_Al*etime(i)))*exp(-lambda_Al*tb);
    line(NBe_p(:,i),NAl_p(:,i)./NBe_p(:,i),'color','k','linestyle',':');    
end