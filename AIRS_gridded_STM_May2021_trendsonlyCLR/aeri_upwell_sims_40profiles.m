%% see aeri_upwell_sims_rtp.m

addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools
addpath /asl/matlib/science/
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations

[h5,ha,p5,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm.rtp');
[p5.salti,p5.landfrac] = usgs_deg10_dem(p5.rlat,p5.rlon);

p5.upwell = 2 * ones(size(p5.upwell));
p5.scanang = 0 * ones(size(p5.scanang));
p5.satzen = 0 * ones(size(p5.scanang));
p5.solzen = 150 * ones(size(p5.scanang));

%%%%%%%%%%%%%%%%%%%%%%%%%

iProf = 10; %% midlats
iProf = 02; %% polar
iProf = 20; %% default, tropics

[hnew,pnew] = replicate_rtp_headprof(h5,p5,iProf,5);  %% raw, increase WV, CO2, T, all three
pnew.upwell = 2 * ones(size(pnew.stemp));
pnew.scanang = 0 * pnew.scanang;
pnew.satzen  = 0 * pnew.satzen;
pnew.solzen  = 150 * ones(size(pnew.satzen));

pnew.gas_1(:,2) = pnew.gas_1(:,2) * 1.001;
pnew.gas_2(:,3) = pnew.gas_2(:,3) * (1 + 2.2/380);
pnew.ptemp(:,4) = pnew.ptemp(:,4) + 0.02;

pnew.gas_1(:,5) = pnew.gas_1(:,5) * 1.001;
pnew.gas_2(:,5) = pnew.gas_2(:,5) * (1 + 2.2/380);
pnew.ptemp(:,5) = pnew.ptemp(:,5) + 0.02;

rtpwrite('upwell.op.rtp',hnew,ha,pnew,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%

%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations

p5pertWV   = p5;   p5pertWV.gas_1   = p5pertWV.gas_1 * (1 + 0.001);
p5pertCO2  = p5;   p5pertCO2.gas_2  = p5pertCO2.gas_2 * (1 + 2.2/385);
p5pertT    = p5;   p5pertT.ptemp    = p5pertT.ptemp + 0.02;
p5pertALL  = p5;   p5pertALL.gas_1 = p5pertALL.gas_1 * (1 + 0.001);
                   p5pertALL.gas_2 = p5pertALL.gas_2 * (1 + 2.2/385);
                   p5pertALL.ptemp = p5pertALL.ptemp + 0.02;

p5pertX = p5pertWV;
[h5,p5pertX] = cat_rtp(h5,p5pertX,h5,p5pertCO2);
[h5,p5pertX] = cat_rtp(h5,p5pertX,h5,p5pertT);
[h5,p5pertX] = cat_rtp(h5,p5pertX,h5,p5pertALL);

p5pertT_TS = p5pertT;     p5pertT_TS.stemp = p5pertT_TS.stemp + 0.02;
p5pertALL_TS = p5pertALL; p5pertALL_TS.stemp = p5pertALL_TS.stemp + 0.02;
p5pertX_TS = p5pertWV;
[h5,p5pertX_TS] = cat_rtp(h5,p5pertX_TS,h5,p5pertCO2);
[h5,p5pertX_TS] = cat_rtp(h5,p5pertX_TS,h5,p5pertT_TS);
[h5,p5pertX_TS] = cat_rtp(h5,p5pertX_TS,h5,p5pertALL_TS);

%rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm_uplook.rtp',h5,ha,p5,pa);
%rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm_uplook_pert.rtp',h5,ha,p5pertX,pa);

mmw = mmwater_rtp(h5,p5);
p5.plays = plevs2plays(p5.plevs);

for ii = 1 : 40
  nlay = p5.nlevs(ii) - 1;
  play = p5.plays(1:nlay,ii);
  tlay = p5.ptemp(1:nlay,ii);
  t2m(ii) = interp1(log(play),tlay,log(p5.spres(ii)),[],'extrap');
  tlowestlay(ii) = tlay(end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the KCARTA DATA
%% get the KCARTA DATA
%% get the KCARTA DATA

%surf_prop_comment: 'stemp mmw co2ppm o3du ch4ppm emiss(900cm-1)  rlat rlon'

iCKD = 1;
iCKD = 25;
iCKD = 32;

fprintf(1,'iCKD = %2i \n',iCKD)

RRTMbands = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];

if ~exist('ilrX')
  disp('need to load 40 x 4 perturbation files , "+" indicates 10, "o" indicates 100')
  dir0 = ['/umbc/xfs2/strow/asl/s1/sergio/home/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/40profiles_uplook_flux_jacs/CKD' num2str(iCKD) '/PERTURBATIONS/'];
  for ii = 1 : 160
    if mod(ii,100) == 0
      fprintf(1,'o');
    elseif mod(ii,10) == 0
      fprintf(1,'+');
    else
      fprintf(1,'.');
    end
    fin = [dir0 '/allkc5bands_prof' num2str(ii) '.mat'];
    junk = load(fin);
    ilrX(ii) = trapz(junk.wall,junk.fluxall(:,1))/1000;  %% to to W/m2
    for bbb = 1 : length(RRTMbands)-1
      boo = find(junk.wall >= RRTMbands(bbb) & junk.wall < RRTMbands(bbb+1));
      ilrXbands(ii,bbb) = trapz(junk.wall(boo),junk.fluxall(boo,1))/1000;  %% to to W/m2
    end
  end
  fprintf(1,'\n');
end

if ~exist('surfprof')
  disp('need to load 40 files , "+" indicates 10')
  dir0 = ['/umbc/xfs2/strow/asl/s1/sergio/home/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/40profiles_uplook_flux_jacs/CKD' num2str(iCKD) '/CALC0_JACOBIANS'];
  for ii = 1 : 40
    if mod(ii,10) == 0
      fprintf(1,'+');
    else
      fprintf(1,'.');
    end
    fin = [dir0 '/allkc5bands_prof' num2str(ii) '.mat'];
    junk = load(fin,'surf_properties');
    surfprof(ii,:) = junk.surf_properties;
  
    junk = load(fin);
    rad0(ii) = trapz(junk.wall,junk.dall);
    ilr0(ii) = trapz(junk.wall,junk.fluxall(:,1))/1000;  %% to to W/m2

    for bbb = 1 : length(RRTMbands)-1
      boo = find(junk.wall >= RRTMbands(bbb) & junk.wall < RRTMbands(bbb+1));
      ilr0bands(ii,bbb) = trapz(junk.wall(boo),junk.fluxall(boo,1))/1000;  %% to to W/m2
    end

    for jjj = 1 : 6
      radpert(ii,jjj) = trapz(junk.wall,junk.jall(:,jjj));    
    end
  end
  fprintf(1,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat = surfprof(:,7);

figure(1); clf

flux0    = pi*rad0'/1000;      %% go from mW/m2/sr-1 to W/m2
fluxpert = pi*radpert/1000;    %% go from mW/m2/sr-1 to W/m2/sr-1

tX2m = surfprof(:,1);
tX2m = t2m';
tX2m = tlowestlay';

plot(lat,ilr0./flux0','r'); hl = legend('\pi rad / ILR','location','best');
plot(lat,flux0,'b',lat,ilr0,'r'); hl = legend('\pi rad','ILR','location','best');

sb = 5.67e-8;   %% stefan botlznann
emiss = flux0 ./ (sb * tX2m.^4);
e2m = (emiss/1.24).^7 .* tX2m;    %% mb

plot(lat,emiss,'b.',lat,1.24*(e2m./tX2m).^(1/7))
plot(lat,flux0,'b.-',lat,sb * emiss .* tX2m.^4)

%% for finite diff jacs, dQ = 0.001, dT = 0.01;
dQ = 0.001;
dT = 0.01;
figure(1); clf; 
plot(1:40,(fluxpert(:,1)-flux0),1:40,fluxpert(:,2:4)-flux0); hl = legend('WV','CO2','O3','CH4','location','best','fontsize',10);

iWV = 4;
iWV = 1;

%iWV

if iWV == 4
  fluxpert(:,1) = fluxpert(:,1) * 4;
  qjac(:,1) = (fluxpert(:,1)-4*flux0)/dQ; %% remember WV has 4 components
elseif iWV == 1
  fluxpert(:,1) = fluxpert(:,1) * 1;
  qjac(:,1) = 4*(fluxpert(:,1)-flux0)/dQ;   %% ah, I fixed KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/cluster_loop_allkcartabands.m
end

qjac(:,2:4) = (fluxpert(:,2:4)-flux0)/dQ;
tjac = (fluxpert(:,5:6)-flux0)/dT;

tpert   = tjac(:,1)*0.02;    %% our paper shows <dST> = 0.02 K/yr
wvpert  = qjac(:,1)*0.001;   %% therefore Eqn 3 of paper (how fracWV changes for const RH) = Lv/Rd dT/T2 = 5422 * 0.02/300/300 = 0.0012
co2pert = qjac(:,2)*2.2/385;
o3pert  = qjac(:,3)*0.001;
ch4pert = qjac(:,4)*6/1800;

figure(1); clf
sumpert = wvpert + co2pert + o3pert + ch4pert + tpert;
plot(lat,wvpert,'b',lat,co2pert,'g',lat,o3pert,lat,ch4pert,lat,tpert,'r',lat,sumpert,'kx-','linewidth',2); 
plotaxis2;
hl = legend('WV','CO2','O3','CH4','T','total','location','best','fontsize',10);
ylabel('Flux change W/m2'); title('for realistic annual change'); xlabel('Latitude')

figure(2); clf
plot(lat,wvpert./sumpert,'b',lat,co2pert./sumpert,'g',lat,o3pert./sumpert,lat,ch4pert./sumpert,lat,tpert./sumpert,'r','linewidth',2); 
plotaxis2;
hl = legend('WV','CO2','O3','CH4','T','location','best','fontsize',10);; xlim([-1 +1]*80);
ylabel('Fractional Contribution \newline to Flux change'); title('for realistic annual change'); xlabel('Latitude');

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf
%% now estimate zonal dependence for CO2
boo = polyfit(lat,co2pert./sumpert,4); sergio_dCO2 = polyval(boo,lat);
boo = polyfit(abs(lat),co2pert./sumpert,1); sergio_dCO2 = polyval(boo,abs(lat));
sergio_dCO2 = 0.1 + 0.2 * sin(pi/180*lat).^2;

%% from Brutsaert d(ILR)/ILR = 27/7 * dT/T + 1/7 dWVfrac/WVfrac
brut_dT  = (27/7 * 0.02./tX2m) .* (sb * emiss .*  tX2m.^4);
brut_dWV = (1/7 * 0.001) .* (sb * emiss .* tX2m.^4);   %% as axel kleidon points out.de/e = 0.07 per kelvin from Clauisus Clapeyron
plot(lat,tpert,'r',lat,brut_dT,'r--',lat,wvpert,'b',lat,brut_dWV,'b--','linewidth',2)
  plotaxis2; ylabel('Flux change W/m2'); xlabel('Latitude')
  hl = legend('tpert actual','tpert Brutsaert','wvpert actual','wvpert Brutsaert','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf
total_dStuff = brut_dT + brut_dWV;
sumpert = wvpert + tpert;
plot(lat,wvpert./sumpert,'b',lat,tpert./sumpert,'r',...
     lat,brut_dWV./total_dStuff,'b--',lat,brut_dT./total_dStuff,'r--','linewidth',2); 
plotaxis2;
hl = legend('WV','T','location','best','fontsize',10);; xlim([-1 +1]*80);
ylabel('Fractional Contribution \newline to Flux change'); title('for realistic annual change'); xlabel('Latitude');

sumpert = wvpert + co2pert + o3pert + ch4pert + tpert;
plot(lat,wvpert./sumpert,'b',lat,co2pert./sumpert,'g',lat,o3pert./sumpert,lat,ch4pert./sumpert,lat,tpert./sumpert,'r',...
     lat,brut_dWV./total_dStuff,'b--',lat,brut_dT./total_dStuff,'r--',lat,sergio_dCO2,'g--','linewidth',2); 
plotaxis2;
hl = legend('WV','CO2','O3','CH4','T','location','best','fontsize',10);; xlim([-1 +1]*80);
ylabel('Fractional Contribution \newline to Flux change'); title('for realistic annual change'); xlabel('Latitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); clf;; 

kcartaWV     = ilrX((1:40)+0*40)-ilr0;
kcartaCO2    = ilrX((1:40)+1*40)-ilr0;
kcartaT      = ilrX((1:40)+2*40)-ilr0;
kcartatotal  = ilrX((1:40)+3*40)-ilr0;
kcartatotal2 = kcartaWV + kcartaCO2 + kcartaT;

plot(lat,ilrX((1:40)+0*40)-ilr0,'b',lat,ilrX((1:40)+1*40)-ilr0,'g',lat,ilrX((1:40)+2*40)-ilr0,'r',lat,ilrX((1:40)+3*40)-ilr0,'kx-','linewidth',2)
plot(lat,kcartaWV,'b',lat,kcartaCO2,'g',lat,kcartaT,'r',lat,kcartatotal,'kx-','linewidth',2)
ylim([0 0.2])
title(['Actual KCARTA Flux Changes!!! CKD ' num2str(iCKD)])
plotaxis2; hl = legend('WV','CO2','T','total');
ylabel('Flux change W/m2'); xlabel('Latitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad_documentation.pdf
% The input NetCDF file contains numerous floating-point variables
% listed in Table 2.1. The dimensions are shown in the order that they
% are listed by the ncdump utility, with the first dimension varying
% slowest in the file. Note that this is opposite to the internal
% Fortran ordering. Most variables are stored as a function of column
% and level (dimensions named col and level in Table 2.1, although the
% actual dimension names are ignored by ecrad). The half level dimension
% corresponds to the mid-points of the levels, plus the
% top-of-atmosphere and surface, and so must be one more than level. The
% level interface dimension excludes the top-of-atmosphere and surface
% so must be one less than level. The optional sw albedo band and lw
% emiss band dimensions allow for shortwave albedo and longwave
% emissivity to be specified in user-defined spectral intervals, but in
% the offline code these are ignored: the first element along these
% dimensions will be used for the entire shortwave and longwave
% spectrum.

%% see /umbc/xfs2/strow/asl/s1/sergio/home/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/wrapper_run_ecRad_rtp_loop_over_profiles.m
if ~exist('ecRad')
  ecRadprofX = superdriver_run_ecRad_rtp_loop_over_profiles(hnew,pnew,-1,-1);
  ecRad0     = superdriver_run_ecRad_rtp_loop_over_profiles(h5,p5,-1,-1);
  %ecRad      = superdriver_run_ecRad_rtp_loop_over_profiles(h5,p5pertX,-1,-1);
  ecRad      = superdriver_run_ecRad_rtp_loop_over_profiles(h5,p5pertX_TS,-1,-1);
end

comment = 'see aeri_upwell_sims_40profles.m';
saver = ['save ilrcalcs_CKD' num2str(iCKD,'%02i') '.mat iCKD ilr0 ilr0bands ilrX ilrXbands ecRad ecRad0 RRTMbands comment']; 
% eval(saver)

ecRadWV     = ecRad.clr((1:40)+0*40)-ecRad0.clr;
ecRadCO2    = ecRad.clr((1:40)+1*40)-ecRad0.clr;
ecRadT      = ecRad.clr((1:40)+2*40)-ecRad0.clr;
ecRadtotal  = ecRad.clr((1:40)+3*40)-ecRad0.clr;
ecRadtotal2 = ecRadWV + ecRadCO2 + ecRadT;

figure(6);
plot(lat,ecRadWV,'b',lat,ecRadCO2,'g',lat,ecRadT,'r',lat,ecRadtotal,'kx-','linewidth',2)
ylim([0 0.2])
title('Actual ecRad Flux Changes!!!')
plotaxis2; hl = legend('WV','CO2','T','total');
ylabel('Flux change W/m2'); xlabel('Latitude')

figure(7)
plot(lat,ecRadWV./ecRadtotal,'b',lat,ecRadCO2./ecRadtotal,'g',lat,ecRadT./ecRadtotal,'r','linewidth',2)
plotaxis2;
hl = legend('WV','CO2','T','location','best','fontsize',10);; xlim([-1 +1]*80);
ylabel('Fractional Contribution \newline to Flux change'); title('for realistic annual change : ecRad'); xlabel('Latitude');

figure(5); ax5 = axis; figure(6); ax6 = axis;
ym = max(ax5(4),ax6(4)); axx = [-75 +75 0 ym]; figure(5); axis(axx); figure(6); axis(axx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); plot([lat; lat; lat; lat],ilrX./ecRad.clr,'linewidth',2);
figure(8); plot(lat,ilr0./ecRad0.clr,'linewidth',2); 
title('Ratio KCARTA Flux/ecRad')

figure(8); plot([lat; lat; lat; lat],ilrX,[lat; lat; lat; lat],ecRad.clr,'linewidth',2);
figure(8); plot(lat,ilr0,lat,ecRad0.clr,'linewidth',2); 
plotaxis2; hl = legend('KCARTA','ecRad'); ylabel('W/m2'); title('Raw Flux')

figure(8); plot([lat; lat; lat; lat],ilrX-[ilr0 ilr0 ilr0 ilr0],'b.',[lat; lat; lat; lat],ecRad.clr-[ecRad0.clr ecRad0.clr ecRad0.clr ecRad0.clr],'r.','linewidth',2);
plotaxis2; title('Perturbations (b) kcarta (r) ecRad')

figure(8); plot(lat,kcartaWV,'b',lat,kcartaCO2,'g',lat,kcartaT,'r',lat,kcartatotal,'k',...
                lat,wvpert,'b:',lat,co2pert,'g:',lat,tpert,'r:',lat,sumpert,'k:',...
                lat,ecRadWV,'b--',lat,ecRadCO2,'g--',lat,ecRadT,'r--',lat,ecRadtotal,'k--','linewidth',2);
plotaxis2; title('Perturbations (solid) kcarta flux \newline (dotted) kcarta jac (dashed) ecRad')
hl = legend('WV','CO2','T','total','location','best','fontsize',10);; xlim([-1 +1]*80);
ylim([0 0.2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
pcolor(p5.rlat,p5.plays,p5.ptemp); caxis([200 300]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap jet; shading interp; ylim([10 1000])
xlabel('Latitude'); ylabel('Pressure (mb)'); colorbar; title('T(p,lat) in K')

figure(10)
ppmv1 = layers2ppmv(h5,p5,1:40,1);
pcolor(p5.rlat,p5.plays(1:97,:),layers2ppmv(h5,p5,1:40,1)); 
caxis([0 3]*1e4); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap jet; shading interp; ylim([100 1000])
xlabel('Latitude'); ylabel('Pressure (mb)'); colorbar; title('WV(p,lat) in g/g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
plot(lat,kcartaWV,'b',lat,kcartaCO2,'g',lat,kcartaT,'r',lat,kcartatotal,'k',lat,wvpert,'b--',lat,co2pert,'g--',lat,tpert,'r--',lat,sumpert,'k--','linewidth',2);
plotaxis2;
hl = legend('WV','CO2','T','total','location','best','fontsize',10);
ylabel('Flux change W/m2'); title('kcarta flux change \newline (a) from flux (dashed) from jacobians')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11);
pcolor(meanvaluebin(RRTMbands),1:160,ilrXbands - [ilr0bands; ilr0bands; ilr0bands; ilr0bands]); colorbar; shading flat;; colormap jet; 
title('KCARTA perturbations'); xlabel('Wavenumber cm-1'); ylabel('WV(1-40),CO2(41-80) \newline T(81-120),All(121-160)')
caxis([0 0.04])

figure(12);
pcolor(meanvaluebin(RRTMbands),1:160,ecRad.bands.clr' - [ecRad0.bands.clr ecRad0.bands.clr ecRad0.bands.clr ecRad0.bands.clr]'); colorbar; shading flat;; colormap jet; 
title('ecRad perturbations'); xlabel('Wavenumber cm-1'); ylabel('WV(1-40),CO2(41-80) \newline T(81-120),All(121-160)')
caxis([0 0.04])

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11);
ind = (1:40)+0*40;
pcolor(meanvaluebin(RRTMbands),1:40,ilrXbands(ind,:) - ilr0bands); colorbar; shading flat;; colormap jet; 
title('KCARTA perturbations'); xlabel('Wavenumber cm-1'); ylabel('WV flux change');
caxis([0 0.04])

figure(12);
pcolor(meanvaluebin(RRTMbands),1:40,ecRad.bands.clr(:,ind)' - ecRad0.bands.clr'); colorbar; shading flat;; colormap jet; 
title('ecRad perturbations'); xlabel('Wavenumber cm-1'); ylabel('WV flux change');
caxis([0 0.04])

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11);
ind = (1:40)+2*40;
pcolor(meanvaluebin(RRTMbands),1:40,ilrXbands(ind,:) - ilr0bands); colorbar; shading flat;; colormap jet; 
title('KCARTA perturbations'); xlabel('Wavenumber cm-1'); ylabel('T flux change');
caxis([0 0.02])

figure(12);
pcolor(meanvaluebin(RRTMbands),1:40,ecRad.bands.clr(:,ind)' - ecRad0.bands.clr'); colorbar; shading flat;; colormap jet; 
title('ecRad perturbations'); xlabel('Wavenumber cm-1'); ylabel('T flux change');
caxis([0 0.02])

