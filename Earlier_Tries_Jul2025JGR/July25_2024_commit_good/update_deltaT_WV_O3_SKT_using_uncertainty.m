if ~exist('results0')
  %% this is the first time code is called, so we can do plots Figs 1,2,3,4,5,6
  results0 = results; 
  deltaT0  = deltaT;
  fracWV0  = fracWV;
  fracO30  = fracO3;
  era50    = era5;

  airsL30     = airsL3;
  cmip60      = cmip6;
  merra20     = merra2;
  climcapsL30 = climcapsL3;
  amip60      = amip6;

end

figure(1); clf
dn = 0:0.01:10; semilogy(dn,histc(resultsunc(:,6)./abs(results(:,6)),dn),'k',dn,histc(deltaTunc(:)./abs(deltaT(:)),dn),'r',dn,histc(fracWVunc(:)./abs(fracWV(:)),dn),'b',dn,histc(fracO3unc(:)./abs(fracO3(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('BEFORE UMBC'); ylabel('\sigma_{unc}/\mu','interpreter','tex')

figure(2); clf
dn = 0:0.01:10; semilogy(dn,histc(era5.trend_stemp_err(:)./abs(era5.trend_stemp(:)),dn),'k',dn,histc(era5.trend_ptemp_err(:)./abs(era5.trend_ptemp(:)),dn),'r',...
                         dn,histc(era5.trend_gas_1_err(:)./abs(era5.trend_gas_1(:)),dn),'b',dn,histc(era5.trend_gas_3_err(:)./abs(era5.trend_gas_3(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('BEFORE ERA5'); ylabel('\sigma_{unc}/\mu','interpreter','tex')

figure(1); clf
dn = 0:0.01:10; semilogy(dn,histc(abs(results(:,6)./results0(:,6)),dn),'k',dn,histc(abs(deltaT(:)./deltaT0(:)),dn),'r',dn,histc(abs(fracWV(:)./fracWV0(:)),dn),'b',dn,histc(abs(fracO3(:)./fracO30(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('BEFORE UMBC'); ylabel('\mu/\mu 0','interpreter','tex')

figure(2); clf
dn = 0:0.01:10; semilogy(dn,histc(abs(era5.trend_stemp(:)./era50.trend_stemp(:)),dn),'k',dn,histc(abs(era5.trend_ptemp(:)./era50.trend_ptemp(:)),dn),'r',...
                         dn,histc(abs(era5.trend_gas_1(:)./era50.trend_gas_1(:)),dn),'b',dn,histc(abs(era5.trend_gas_3(:)./era50.trend_gas_3(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('BEFORE ERA5'); ylabel('\mu/\mu 0','interpreter','tex')

maxratio = 2;
maxratio = 0.25;
maxratio = 0.50;
maxratio = 0.75;
maxratio = 1.00;

disp(' ')
fprintf(1,'updating deltaT/deltaWV/deltaO3/deltaSKT according to uncertainty : max delta factor allowed = %8.3f \n',maxratio);
fprintf(1,'updating deltaT/deltaWV/deltaO3/deltaSKT according to uncertainty : max delta factor allowed = %8.3f \n',maxratio);
fprintf(1,'updating deltaT/deltaWV/deltaO3/deltaSKT according to uncertainty : max delta factor allowed = %8.3f \n',maxratio);
disp(' ')

if exist('results0') & exist('era50')
  %% if it don't work, try and try again!!!!
  results    = results0; 
  deltaT     = deltaT0;
  fracWV     = fracWV0;
  fracO3     = fracO30;
  era5       = era50;
  airsL3     = airsL30;
  cmip6      = cmip60;
  merra2     = merra20;
  climcapsL3 = climcapsL30;
  amip6      = amip60;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results0 = results; 
deltaT0  = deltaT;
fracWV0  = fracWV;
fracO30  = fracO3;

figure(3); clf; plot(resultsunc(:,6) ./ results(:,6)); title('BEFORE deltaSTunc/deltaST')
figure(4); clf; subplot(121); pcolor(deltaT); set(gca,'ydir','reverse'); shading flat; colorbar; colormap(usa2); title('BEFORE deltaT'); ;  caxis([-1 +1]*0.25)
                subplot(122); pcolor(abs(deltaTunc./deltaT)); set(gca,'ydir','reverse'); shading flat; colorbar; title('BEFORE unc/delta'); caxis([0 +1])
figure(5); clf; subplot(121); pcolor(fracWV); set(gca,'ydir','reverse'); shading flat; colorbar; colormap(usa2); title('BEFORE fracWV');    caxis([-1 +1]*0.025)
                subplot(122); pcolor(abs(fracWVunc./fracWV)); set(gca,'ydir','reverse'); shading flat; colorbar; title('BEFORE unc/delta'); caxis([0 +1])
figure(6); clf; subplot(121); pcolor(fracO3); set(gca,'ydir','reverse'); shading flat; colorbar; colormap(usa2); title('BEFORE fracO3');    caxis([-1 +1]*0.015)
                subplot(122); pcolor(abs(fracO3unc./fracO3)); set(gca,'ydir','reverse'); shading flat; colorbar; title('BEFORE unc/delta'); caxis([0 +1])

ssgn = sign(results(:,6));
ssgn = ones(size(ssgn));
% bad = resultsunc(:,6)./abs(results(:,6)); % bad = find(bad > maxratio); resultsunc(bad,6) = abs(results(bad,6))*maxratio;
% bad = resultsunc(:,6)./abs(results(:,6)); % bad = find(bad > maxratio); resultsunc(bad,6) = abs(results(bad,6))*maxratio;
%   results(:,6) = results(:,6) + ssgn .* resultsunc(:,6);
junk = abs(results(:,6))*maxratio; bad=find(resultsunc(:,6) > junk);  junk = min(junk,resultsunc(:,6));
  results(:,6) = results(:,6) + ssgn .* junk;
  fprintf(1,'for UMBC   : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad),length(junk),length(bad)/length(junk))

ssgn = sign(deltaT);
ssgn = ones(size(ssgn));
% bad = deltaTunc./abs(deltaT); % bad = find(bad > maxratio); deltaTunc(bad) = abs(deltaT(bad))*maxratio;
% bad = deltaTunc./abs(deltaT); % bad = find(bad > maxratio); deltaTunc(bad) = abs(deltaT(bad))*maxratio;
%   deltaT = deltaT + ssgn.*deltaTunc;
junk = abs(deltaT)*maxratio; bad=find(deltaTunc > junk);  junk = min(junk,deltaTunc);
  deltaT = deltaT + ssgn.*junk;
  fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(fracWV);
ssgn = ones(size(ssgn));
% bad = fracWVunc./abs(fracWV); % bad = find(bad > maxratio); fracWVunc(bad) = abs(fracWV(bad))*maxratio;
% bad = fracWVunc./abs(fracWV); % bad = find(bad > maxratio); fracWVunc(bad) = abs(fracWV(bad))*maxratio;
%   fracWV = fracWV + ssgn.*fracWVunc;
junk = abs(fracWV)*maxratio; bad=find(fracWVunc > junk);  junk = min(junk,fracWVunc);
  fracWV = fracWV + ssgn.*junk;
  fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(fracO3);
ssgn = ones(size(ssgn));
% bad = fracO3unc./abs(fracO3); % bad = find(bad > maxratio); fracO3unc(bad) = abs(fracO3(bad))*maxratio;
% bad = fracO3unc./abs(fracO3); % bad = find(bad > maxratio); fracO3unc(bad) = abs(fracO3(bad))*maxratio;
%   fracO3 = fracO3 + ssgn.*fracO3unc;
junk = abs(fracO3)*maxratio; bad=find(fracO3unc > junk);  junk = min(junk,fracO3unc);
  fracO3 = fracO3 + ssgn.*junk;
  fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

%%%%%%%%%%%%%%%%%%%%%%%%%

era50 = era5;

ssgn = sign(era5.trend_stemp);
ssgn = ones(size(ssgn));
% bad = era5.trend_stemp_err ./ abs(era5.trend_stemp); % bad = find(bad > maxratio); era5.trend_stemp_err(bad) = abs(era5.trend_stemp(bad))*maxratio;
%   era5.trend_stemp = era5.trend_stemp + ssgn.*era5.trend_stemp_err;
junk = abs(era5.trend_stemp)*maxratio; bad=find(era5.trend_stemp_err > junk);  junk = min(junk,era5.trend_stemp_err);
  era5.trend_stemp = era5.trend_stemp + ssgn.*junk;
  fprintf(1,'for ERA5   : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad),length(junk),length(bad)/length(junk))

ssgn = sign(era5.trend_ptemp);
ssgn = ones(size(ssgn));
% bad = era5.trend_ptemp_err ./ abs(era5.trend_ptemp); % bad = find(bad > maxratio); era5.trend_ptemp_err(bad) = abs(era5.trend_ptemp(bad))*maxratio;
%   era5.trend_ptemp = era5.trend_ptemp + ssgn.*era5.trend_ptemp_err;
junk = abs(era5.trend_ptemp)*maxratio; bad=find(era5.trend_ptemp_err > junk);  junk = min(junk,era5.trend_ptemp_err);
  era5.trend_ptemp = era5.trend_ptemp + ssgn.*junk;
  fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(era5.trend_gas_1);
ssgn = ones(size(ssgn));
% bad = era5.trend_gas_1_err ./ abs(era5.trend_gas_1); % bad = find(bad > maxratio); era5.trend_gas_1_err(bad) = abs(era5.trend_gas_1(bad))*maxratio;
%   era5.trend_gas_1 = era5.trend_gas_1 + ssgn.*era5.trend_gas_1_err;
junk = abs(era5.trend_gas_1)*maxratio; bad=find(era5.trend_gas_1_err > junk);  junk = min(junk,era5.trend_gas_1_err);
  era5.trend_gas_1 = era5.trend_gas_1 + ssgn.*junk;
  fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(era5.trend_gas_3);
ssgn = ones(size(ssgn));
% bad = era5.trend_gas_3_err ./ abs(era5.trend_gas_3); % bad = find(bad > maxratio); era5.trend_gas_3_err(bad) = abs(era5.trend_gas_3(bad))*maxratio;
%   era5.trend_gas_3 = era5.trend_gas_3 + ssgn.*era5.trend_gas_3_err;
junk = abs(era5.trend_gas_3)*maxratio; bad=find(era5.trend_gas_3_err > junk);  junk = min(junk,era5.trend_gas_3_err);
  era5.trend_gas_3 = era5.trend_gas_3 + ssgn.*junk;
  fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

%%%%%%%%%%%%%%%%%%%%%%%%%

airsL30 = airsL3;

ssgn = sign(airsL3.thestats64x72.stemprate);
ssgn = ones(size(ssgn));
% bad = airsL3.thestats64x72.stempratestd ./ abs(airsL3.thestats64x72.stemprate); % bad = find(bad > maxratio); airsL3.thestats64x72.stempratestd(bad) = abs(airsL3.thestats64x72.stemprate(bad))*maxratio;
%   airsL3.thestats64x72.stemprate = airsL3.thestats64x72.stemprate + ssgn.*airsL3.thestats64x72.stempratestd;
junk = abs(airsL3.thestats64x72.stemprate)*maxratio; bad=find(airsL3.thestats64x72.stempratestd > junk);  junk = min(junk,airsL3.thestats64x72.stempratestd);
  airsL3.thestats64x72.stemprate = airsL3.thestats64x72.stemprate + ssgn.*junk;
  fprintf(1,'for AIRSL3 : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(airsL3.thestats64x72.ptemprate);
ssgn = ones(size(ssgn));
% bad = airsL3.thestats64x72.ptempratestd ./ abs(airsL3.thestats64x72.ptemprate); % bad = find(bad > maxratio); airsL3.thestats64x72.ptempratestd(bad) = abs(airsL3.thestats64x72.ptemprate(bad))*maxratio;
%   airsL3.thestats64x72.ptemprate = airsL3.thestats64x72.ptemprate + ssgn.*airsL3.thestats64x72.ptempratestd;
junk = abs(airsL3.thestats64x72.ptemprate)*maxratio; bad=find(airsL3.thestats64x72.ptempratestd > junk);  junk = min(junk,airsL3.thestats64x72.ptempratestd);
  airsL3.thestats64x72.ptemprate = airsL3.thestats64x72.ptemprate + ssgn.*junk;
  fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(airsL3.thestats64x72.waterrate);
ssgn = ones(size(ssgn));
% bad = airsL3.thestats64x72.waterratestd ./ abs(airsL3.thestats64x72.waterrate); % bad = find(bad > maxratio); airsL3.thestats64x72.waterratestd(bad) = abs(airsL3.thestats64x72.waterrate(bad))*maxratio;
%   airsL3.thestats64x72.waterrate = airsL3.thestats64x72.waterrate + ssgn.*airsL3.thestats64x72.waterratestd;
junk = abs(airsL3.thestats64x72.waterrate)*maxratio; bad=find(airsL3.thestats64x72.waterratestd > junk);  junk = min(junk,airsL3.thestats64x72.waterratestd);
  airsL3.thestats64x72.waterrate = airsL3.thestats64x72.waterrate + ssgn.*junk;
  fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(airsL3.thestats64x72.ozonerate);
ssgn = ones(size(ssgn));
% bad = airsL3.thestats64x72.ozoneratestd ./ abs(airsL3.thestats64x72.ozonerate); % bad = find(bad > maxratio); airsL3.thestats64x72.ozoneratestd(bad) = abs(airsL3.thestats64x72.ozonerate(bad))*maxratio;
%   airsL3.thestats64x72.ozonerate = airsL3.thestats64x72.ozonerate + ssgn.*airsL3.thestats64x72.ozoneratestd;
junk = abs(airsL3.thestats64x72.ozonerate)*maxratio; bad=find(airsL3.thestats64x72.ozoneratestd > junk);  junk = min(junk,airsL3.thestats64x72.ozoneratestd);
  airsL3.thestats64x72.ozonerate = airsL3.thestats64x72.ozonerate + ssgn.*junk;
  fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

%%%%%%%%%%%%%%%%%%%%%%%%%

cmip60 = cmip6;

ssgn = sign(cmip6.trend_stemp);
ssgn = ones(size(ssgn));
% bad = cmip6.trend_stemp_err ./ abs(cmip6.trend_stemp); % bad = find(bad > maxratio); cmip6.trend_stemp_err(bad) = abs(cmip6.trend_stemp(bad))*maxratio;
%   cmip6.trend_stemp = cmip6.trend_stemp + ssgn.*cmip6.trend_stemp_err;
junk = abs(cmip6.trend_stemp)*maxratio; bad=find(cmip6.trend_stemp_err > junk);  junk = min(junk,cmip6.trend_stemp_err);
  cmip6.trend_stemp = cmip6.trend_stemp + ssgn.*junk;
  fprintf(1,'for CMIP6  : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad),length(junk),length(bad)/length(junk))

ssgn = sign(cmip6.trend_ptemp);
ssgn = ones(size(ssgn));
% bad = cmip6.trend_ptemp_err ./ abs(cmip6.trend_ptemp); % bad = find(bad > maxratio); cmip6.trend_ptemp_err(bad) = abs(cmip6.trend_ptemp(bad))*maxratio;
%   cmip6.trend_ptemp = cmip6.trend_ptemp + ssgn.*cmip6.trend_ptemp_err;
junk = abs(cmip6.trend_ptemp)*maxratio; bad=find(cmip6.trend_ptemp_err > junk);  junk = min(junk,cmip6.trend_ptemp_err);
  cmip6.trend_ptemp = cmip6.trend_ptemp + ssgn.*junk;
  fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(cmip6.trend_gas_1);
ssgn = ones(size(ssgn));
% bad = cmip6.trend_gas_1_err ./ abs(cmip6.trend_gas_1); % bad = find(bad > maxratio); cmip6.trend_gas_1_err(bad) = abs(cmip6.trend_gas_1(bad))*maxratio;
%   cmip6.trend_gas_1 = cmip6.trend_gas_1 + ssgn.*cmip6.trend_gas_1_err;
junk = abs(cmip6.trend_gas_1)*maxratio; bad=find(cmip6.trend_gas_1_err > junk);  junk = min(junk,cmip6.trend_gas_1_err);
  cmip6.trend_gas_1 = cmip6.trend_gas_1 + ssgn.*junk;
  fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

ssgn = sign(cmip6.trend_gas_3);
ssgn = ones(size(ssgn));
% bad = cmip6.trend_gas_3_err ./ abs(cmip6.trend_gas_3); % bad = find(bad > maxratio); cmip6.trend_gas_3_err(bad) = abs(cmip6.trend_gas_3(bad))*maxratio;
%   cmip6.trend_gas_3 = cmip6.trend_gas_3 + ssgn.*cmip6.trend_gas_3_err;
junk = abs(cmip6.trend_gas_3)*maxratio; bad=find(cmip6.trend_gas_3_err > junk);  junk = min(junk,cmip6.trend_gas_3_err);
  cmip6.trend_gas_3 = cmip6.trend_gas_3 + ssgn.*junk;
  fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))
  
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('merra2')

  merra20 = merra2;

  ssgn = sign(merra2.trend_stemp);
  ssgn = ones(size(ssgn));
  % bad = merra2.trend_stemp_err ./ abs(merra2.trend_stemp); % bad = find(bad > maxratio); merra2.trend_stemp_err(bad) = abs(merra2.trend_stemp(bad))*maxratio;
  %   merra2.trend_stemp = merra2.trend_stemp + ssgn.*merra2.trend_stemp_err;
  junk = abs(merra2.trend_stemp)*maxratio; bad=find(merra2.trend_stemp_err > junk);  junk = min(junk,merra2.trend_stemp_err);
    merra2.trend_stemp = merra2.trend_stemp + ssgn.*junk;
    fprintf(1,'for MERRA2 : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad),length(junk),length(bad)/length(junk))

  ssgn = sign(merra2.trend_ptemp);
  ssgn = ones(size(ssgn));
  % bad = merra2.trend_ptemp_err ./ abs(merra2.trend_ptemp); % bad = find(bad > maxratio); merra2.trend_ptemp_err(bad) = abs(merra2.trend_ptemp(bad))*maxratio;
  %   merra2.trend_ptemp = merra2.trend_ptemp + ssgn.*merra2.trend_ptemp_err;
  junk = abs(merra2.trend_ptemp)*maxratio; bad=find(merra2.trend_ptemp_err > junk);  junk = min(junk,merra2.trend_ptemp_err);
    merra2.trend_ptemp = merra2.trend_ptemp + ssgn.*junk;
    fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(merra2.trend_gas_1);
  ssgn = ones(size(ssgn));
  % bad = merra2.trend_gas_1_err ./ abs(merra2.trend_gas_1); % bad = find(bad > maxratio); merra2.trend_gas_1_err(bad) = abs(merra2.trend_gas_1(bad))*maxratio;
  %   merra2.trend_gas_1 = merra2.trend_gas_1 + ssgn.*merra2.trend_gas_1_err;
  junk = abs(merra2.trend_gas_1)*maxratio; bad=find(merra2.trend_gas_1_err > junk);  junk = min(junk,merra2.trend_gas_1_err);
    merra2.trend_gas_1 = merra2.trend_gas_1 + ssgn.*junk;
    fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(merra2.trend_gas_3);
  ssgn = ones(size(ssgn));
  % bad = merra2.trend_gas_3_err ./ abs(merra2.trend_gas_3); % bad = find(bad > maxratio); merra2.trend_gas_3_err(bad) = abs(merra2.trend_gas_3(bad))*maxratio;
  %   merra2.trend_gas_3 = merra2.trend_gas_3 + ssgn.*merra2.trend_gas_3_err;
  junk = abs(merra2.trend_gas_3)*maxratio; bad=find(merra2.trend_gas_3_err > junk);  junk = min(junk,merra2.trend_gas_3_err);
    merra2.trend_gas_3 = merra2.trend_gas_3 + ssgn.*junk;
    fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  climcapsL30 = climcapsL3;

  ssgn = sign(climcapsL3.thestats64x72.stemprate);
  ssgn = ones(size(ssgn));
  % bad = climcapsL3.thestats64x72.stempratestd ./ abs(climcapsL3.thestats64x72.stemprate); % bad = find(bad > maxratio); climcapsL3.thestats64x72.stempratestd(bad) = abs(climcapsL3.thestats64x72.stemprate(bad))*maxratio;
  %   climcapsL3.thestats64x72.stemprate = climcapsL3.thestats64x72.stemprate + ssgn.*climcapsL3.thestats64x72.stempratestd;
  junk = abs(climcapsL3.thestats64x72.stemprate)*maxratio; bad=find(climcapsL3.thestats64x72.stempratestd > junk);  junk = min(junk,climcapsL3.thestats64x72.stempratestd);
    climcapsL3.thestats64x72.stemprate = climcapsL3.thestats64x72.stemprate + ssgn.*junk;
    fprintf(1,'for CLIML3 : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(climcapsL3.thestats64x72.ptemprate);
  ssgn = ones(size(ssgn));
  % bad = climcapsL3.thestats64x72.ptempratestd ./ abs(climcapsL3.thestats64x72.ptemprate); % bad = find(bad > maxratio); climcapsL3.thestats64x72.ptempratestd(bad) = abs(climcapsL3.thestats64x72.ptemprate(bad))*maxratio;
  %   climcapsL3.thestats64x72.ptemprate = climcapsL3.thestats64x72.ptemprate + ssgn.*climcapsL3.thestats64x72.ptempratestd;
  junk = abs(climcapsL3.thestats64x72.ptemprate)*maxratio; bad=find(climcapsL3.thestats64x72.ptempratestd > junk);  junk = min(junk,climcapsL3.thestats64x72.ptempratestd);
    climcapsL3.thestats64x72.ptemprate = climcapsL3.thestats64x72.ptemprate + ssgn.*junk;
    fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(climcapsL3.thestats64x72.waterrate);
  ssgn = ones(size(ssgn));
  % bad = climcapsL3.thestats64x72.waterratestd ./ abs(climcapsL3.thestats64x72.waterrate); % bad = find(bad > maxratio); climcapsL3.thestats64x72.waterratestd(bad) = abs(climcapsL3.thestats64x72.waterrate(bad))*maxratio;
  %   climcapsL3.thestats64x72.waterrate = climcapsL3.thestats64x72.waterrate + ssgn.*climcapsL3.thestats64x72.waterratestd;
  junk = abs(climcapsL3.thestats64x72.waterrate)*maxratio; bad=find(climcapsL3.thestats64x72.waterratestd > junk);  junk = min(junk,climcapsL3.thestats64x72.waterratestd);
    climcapsL3.thestats64x72.waterrate = climcapsL3.thestats64x72.waterrate + ssgn.*junk;
    fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(climcapsL3.thestats64x72.ozonerate);
  ssgn = ones(size(ssgn));
  % bad = climcapsL3.thestats64x72.ozoneratestd ./ abs(climcapsL3.thestats64x72.ozonerate); % bad = find(bad > maxratio); climcapsL3.thestats64x72.ozoneratestd(bad) = abs(climcapsL3.thestats64x72.ozonerate(bad))*maxratio;
  %   climcapsL3.thestats64x72.ozonerate = climcapsL3.thestats64x72.ozonerate + ssgn.*climcapsL3.thestats64x72.ozoneratestd;
  junk = abs(climcapsL3.thestats64x72.ozonerate)*maxratio; bad=find(climcapsL3.thestats64x72.ozoneratestd > junk);  junk = min(junk,climcapsL3.thestats64x72.ozoneratestd);
    climcapsL3.thestats64x72.ozonerate = climcapsL3.thestats64x72.ozonerate + ssgn.*junk;
    fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  amip60 = amip6;

  ssgn = sign(amip6.trend_stemp);
  ssgn = ones(size(ssgn));
  % bad = amip6.trend_stemp_err ./ abs(amip6.trend_stemp); % bad = find(bad > maxratio); amip6.trend_stemp_err(bad) = abs(amip6.trend_stemp(bad))*maxratio;
  %   amip6.trend_stemp = amip6.trend_stemp + ssgn.*amip6.trend_stemp_err;
  junk = abs(amip6.trend_stemp)*maxratio; bad=find(amip6.trend_stemp_err > junk);  junk = min(junk,amip6.trend_ptemp_err);
    amip6.trend_ptemp = amip6.trend_ptemp + ssgn.*junk;
    fprintf(1,'for AMIP6  : stemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad),length(junk),length(bad)/length(junk))

  ssgn = sign(amip6.trend_ptemp);
  ssgn = ones(size(ssgn));
  % bad = amip6.trend_ptemp_err ./ abs(amip6.trend_ptemp); % bad = find(bad > maxratio); amip6.trend_ptemp_err(bad) = abs(amip6.trend_ptemp(bad))*maxratio;
  %   amip6.trend_ptemp = amip6.trend_ptemp + ssgn.*amip6.trend_ptemp_err;
  junk = abs(amip6.trend_ptemp)*maxratio; bad=find(amip6.trend_ptemp_err > junk);  junk = min(junk,amip6.trend_ptemp_err);
    amip6.trend_ptemp = amip6.trend_ptemp + ssgn.*junk;
    fprintf(1,'           : ptemp : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(amip6.trend_gas_1);
  ssgn = ones(size(ssgn));
  % bad = amip6.trend_gas_1_err ./ abs(amip6.trend_gas_1); % bad = find(bad > maxratio); amip6.trend_gas_1_err(bad) = abs(amip6.trend_gas_1(bad))*maxratio;
  %   amip6.trend_gas_1 = amip6.trend_gas_1 + ssgn.*amip6.trend_gas_1_err;
  junk = abs(amip6.trend_gas_1)*maxratio; bad=find(amip6.trend_gas_1_err > junk);  junk = min(junk,amip6.trend_gas_1_err);
    amip6.trend_gas_1 = amip6.trend_gas_1 + ssgn.*junk;
    fprintf(1,'           : water : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))

  ssgn = sign(amip6.trend_gas_3);
  ssgn = ones(size(ssgn));
  % bad = amip6.trend_gas_3_err ./ abs(amip6.trend_gas_3); % bad = find(bad > maxratio); amip6.trend_gas_3_err(bad) = abs(amip6.trend_gas_3(bad))*maxratio;
  %   amip6.trend_gas_3 = amip6.trend_gas_3 + ssgn.*amip6.trend_gas_3_err;
  junk = abs(amip6.trend_gas_3)*maxratio; bad=find(amip6.trend_gas_3_err > junk);  junk = min(junk,amip6.trend_gas_3_err);
    amip6.trend_gas_3 = amip6.trend_gas_3 + ssgn.*junk;
    fprintf(1,'           : ozone : found %10i of %10i = %8.6f fraction were bad \n',length(bad(:)),length(junk(:)),length(bad(:))/length(junk(:)))
  %%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7); clf
dn = 0:0.01:10; semilogy(dn,histc(resultsunc(:,6)./abs(results(:,6)),dn),'k',dn,histc(deltaTunc(:)./abs(deltaT(:)),dn),'r',dn,histc(fracWVunc(:)./abs(fracWV(:)),dn),'b',dn,histc(fracO3unc(:)./abs(fracO3(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('AFTER UMBC'); ylabel('\sigma_{unc}/\mu','interpreter','tex')

figure(8); clf
dn = 0:0.01:10; semilogy(dn,histc(era5.trend_stemp_err(:)./abs(era5.trend_stemp(:)),dn),'k',dn,histc(era5.trend_ptemp_err(:)./abs(era5.trend_ptemp(:)),dn),'r',...
                         dn,histc(era5.trend_gas_1_err(:)./abs(era5.trend_gas_1(:)),dn),'b',dn,histc(era5.trend_gas_3_err(:)./abs(era5.trend_gas_3(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('AFTER ERA5'); ylabel('\sigma_{unc}/\mu','interpreter','tex')

figure(7); clf
dn = 0:0.01:10; semilogy(dn,histc(abs(results(:,6)./results0(:,6)),dn),'k',dn,histc(abs(deltaT(:)./deltaT0(:)),dn),'r',dn,histc(abs(fracWV(:)./fracWV0(:)),dn),'b',dn,histc(abs(fracO3(:)./fracO30(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('AFTER UMBC'); ylabel('\mu/\mu 0','interpreter','tex')

figure(8); clf
dn = 0:0.01:10; semilogy(dn,histc(abs(era5.trend_stemp(:)./era50.trend_stemp(:)),dn),'k',dn,histc(abs(era5.trend_ptemp(:)./era50.trend_ptemp(:)),dn),'r',...
                         dn,histc(abs(era5.trend_gas_1(:)./era50.trend_gas_1(:)),dn),'b',dn,histc(abs(era5.trend_gas_3(:)./era50.trend_gas_3(:)),dn),'g'); grid;
hl = legend('GHG+Stemp','T(z)','WV(z)','O3(z)','location','best'); title('AFTER ERA5'); ylabel('\mu/\mu 0','interpreter','tex')

figure(9); clf; plot(resultsunc(:,6) ./ results(:,6)); title('AFTER deltaSTunc/deltaST')
figure(10); clf; subplot(121); pcolor(deltaT); set(gca,'ydir','reverse'); shading flat; colorbar; colormap(usa2); title('AFTER deltaT'); ;  caxis([-1 +1]*0.25)
                subplot(122); pcolor(abs(deltaTunc./deltaT)); set(gca,'ydir','reverse'); shading flat; colorbar; title('AFTER unc/delta'); caxis([0 +1])
figure(11); clf; subplot(121); pcolor(fracWV); set(gca,'ydir','reverse'); shading flat; colorbar; colormap(usa2); title('AFTER fracWV');    caxis([-1 +1]*0.025)
                subplot(122); pcolor(abs(fracWVunc./fracWV)); set(gca,'ydir','reverse'); shading flat; colorbar; title('AFTER unc/delta'); caxis([0 +1])
figure(12); clf; subplot(121); pcolor(fracO3); set(gca,'ydir','reverse'); shading flat; colorbar; colormap(usa2); title('AFTER fracO3');    caxis([-1 +1]*0.015)
                subplot(122); pcolor(abs(fracO3unc./fracO3)); set(gca,'ydir','reverse'); shading flat; colorbar; title('AFTER unc/delta'); caxis([0 +1])
pause(0.1)
