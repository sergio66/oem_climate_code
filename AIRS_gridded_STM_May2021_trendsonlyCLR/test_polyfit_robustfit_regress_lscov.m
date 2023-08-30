iNumYears = 20;
if ~exist('p')
  read_fileMean17years
  h = hMean17years;
  p = pMean17years;
end

planet  = 1 : 4608;
tropics = find(abs(p.rlat) <= 30);
poles   = find(abs(p.rlat) > 60);
midlats = setdiff(planet,union(tropics,poles));

figure(1); clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('umbc_spectral_olr')
  feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];
  loader = ['load ' feedbacknameUMBC];
  eval(loader)
end

plot(p.stemp,umbc_spectral_olr.olr0_ecRad.clr,'.')

xUMBC = p.stemp;
yUMBC = umbc_spectral_olr.olr0_ecRad.clr;

xUMBC = results(:,6)';
yUMBC = umbc_spectral_olr.olr0_ecRad.clr - umbc_spectral_olr.wv_ecRad.clr;

[nn,nx,ny,nmean_tropics_UMBC,nstd_tropics_UMBC] = myhist2d(xUMBC(tropics),yUMBC(tropics),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_tropics_UMBC,nstd_tropics_UMBC)

[nn,nx,ny,nmean_midlats_UMBC,nstd_midlats_UMBC] = myhist2d(xUMBC(midlats),yUMBC(midlats),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_midlats_UMBC,nstd_midlats_UMBC)

[nn,nx,ny,nmean_poles_UMBC,nstd_poles_UMBC] = myhist2d(xUMBC(poles),yUMBC(poles),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_poles_UMBC,nstd_poles_UMBC)

[nn,nx,ny,nmean_planet_UMBC,nstd_planet_UMBC] = myhist2d(xUMBC(planet),yUMBC(planet),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_planet_UMBC,nstd_planet_UMBC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('era5_spectral_olr')
  strMODELS = 'AIRSL3_ERA5_CMIP6';
  feedbacknameNWP_ERA5 = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELS '_numyears_' num2str(iNumYears,'%02d') '.mat'];
  loader = ['load ' feedbacknameNWP_ERA5];
  eval(loader)
end

xERA5 = stemptrend.era5;
yERA5 = era5_spectral_olr.olr0_ecRad.clr - era5_spectral_olr.wv_ecRad.clr;

[nn,nx,ny,nmean_tropics_ERA5,nstd_tropics_ERA5] = myhist2d(xERA5(tropics),yERA5(tropics),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_tropics_ERA5,nstd_tropics_ERA5)

[nn,nx,ny,nmean_midlats_ERA5,nstd_midlats_ERA5] = myhist2d(xERA5(midlats),yERA5(midlats),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_midlats_ERA5,nstd_midlats_ERA5)

[nn,nx,ny,nmean_poles_ERA5,nstd_poles_ERA5] = myhist2d(xERA5(poles),yERA5(poles),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_poles_ERA5,nstd_poles_ERA5)

[nn,nx,ny,nmean_planet_ERA5,nstd_planet_ERA5] = myhist2d(xERA5(planet),yERA5(planet),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmean_planet_ERA5,nstd_planet_ERA5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf; 
errorbar((-1:0.01:+1)-0.001,nmean_tropics_UMBC,nstd_tropics_UMBC,'color','b','linewidth',2); hold on
errorbar((-1:0.01:+1)+0.001,nmean_tropics_ERA5,nstd_tropics_ERA5,'color','r','linewidth',2); hold off

xlim([-1 +1]*0.15); plotaxis2; 
hl = legend('UMBC','ERA5','location','best'); 

xlabel('\delta SKT (K)'); ylabel('\delta OLR W/m2');
title('Tropical WV perturbations')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf; 
errorbar((-1:0.01:+1)-0.001,nmean_midlats_UMBC,nstd_midlats_UMBC,'color','b','linewidth',2); hold on
errorbar((-1:0.01:+1)+0.001,nmean_midlats_ERA5,nstd_midlats_ERA5,'color','r','linewidth',2); hold off

xlim([-1 +1]*0.15); plotaxis2; 
hl = legend('UMBC','ERA5','location','best'); 

xlabel('\delta SKT (K)'); ylabel('\delta OLR W/m2');
title('Midlats WV perturbations')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf; 
errorbar((-1:0.01:+1)-0.001,nmean_poles_UMBC,nstd_poles_UMBC,'color','b','linewidth',2); hold on
errorbar((-1:0.01:+1)+0.001,nmean_poles_ERA5,nstd_poles_ERA5,'color','r','linewidth',2); hold off

xlim([-1 +1]*0.15); plotaxis2; 
hl = legend('UMBC','ERA5','location','best'); 

xlabel('\delta SKT (K)'); ylabel('\delta OLR W/m2');
title('Poles WV perturbations')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf; 
errorbar((-1:0.01:+1)-0.001,nmean_planet_UMBC,nstd_planet_UMBC,'color','b','linewidth',2); hold on
errorbar((-1:0.01:+1)+0.001,nmean_planet_ERA5,nstd_planet_ERA5,'color','r','linewidth',2); hold off

xlim([-1 +1]*0.15); plotaxis2; 
hl = legend('UMBC','ERA5','location','best'); 

xlabel('\delta SKT (K)'); ylabel('\delta OLR W/m2');
title('Global WV perturbations')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zone = planet;
zone = tropics; 
zone = poles; 

x = xERA5(zone);
y = yERA5(zone);

[P0]    = polyfit(x,y,1);
[P1 B1] = robustfit(x,y);

modelFun = @(b,x) b(1)*x + b(2);
nlm0 = fitnlm(x,y,modelFun,P0);
nlm1 = fitnlm(x,y,modelFun,P0,'Weight',cos(p.rlat(zone) * pi/180));
nlm2 = fitnlm(x,y,modelFun,P0,'Weight',1./cos(p.rlat(zone) * pi/180));

lm0 = fitlm(x,y);
lm1 = fitlm(x,y,'Weight',cos(p.rlat(zone) * pi/180));
lm2 = fitlm(x,y,'Weight',1./cos(p.rlat(zone) * pi/180));

x1 = ones(size(x));
X  = [x1; x];
[rb,rbint,rr,rrint,rstats] = regress(y',X');

%%%%%%%%%%
fprintf(1,'polyfit      : delta OLR = %8.6f delta SKT + %8.6f \n',P0(1),P0(2))
fprintf(1,'robustfitfit : delta OLR = %8.6f delta SKT + %8.6f \n',P1(2),P1(1))
fprintf(1,'regress      : delta OLR = %8.6f delta SKT + %8.6f \n',rb(2),rb(1));
[nn,nx,ny,nmean,nstd] = myhist2d(x,y,-1:0.01:+1,-1:0.01:+1);
figure(5); clf; plot(x,y,'b.',x,polyval(P0,x),'ro',x,P0(1)*x+P0(2),'g'); hold on
errorbar(-1:0.01:+1,nmean,nstd,'linewidth',2,'color','k'); hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


