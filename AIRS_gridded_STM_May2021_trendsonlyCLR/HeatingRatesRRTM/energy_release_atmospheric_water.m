addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools

[h,ha,p,pa] = rtpread('pertO3.rtp');

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h,p,1:length(p.stemp),1); 
if ~isfield(p,'plays')
  p.plays = plevs2plays(p.plevs);
end
y = ppmv2gg(p.plays(1:97,:),p.ptemp(1:97,:),ppmvLAY,18);

figure(1)
pcolor(y); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar

pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),y); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('WV g/g');

pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),p.gas_1(1:97,:)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('WV molecules/cm2');

mass_w = p.gas_1(1:97,:)*1e4*18*1.67e-27; 
mass_a = mass_w ./ y;

pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),mass_w);
set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('WV kg/m2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
%% we know dT/dt ~ 0.025 K/yr to 100 mb, and -0.025 K/yr above that
dq = zeros(97,40);
boo = find(p.plays(1:97,20) >= 100,1); 
atm_dq_from_dTzdt_per_year(1:boo-1,:)  = mass_a(1:boo-1,:) * 1003.5 * -0.025;
atm_dq_from_dTzdt_per_year(boo:97,:)   = mass_a(boo:97,:)  * 1003.5 * +0.025;
pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),atm_dq_from_dTzdt_per_year);
set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('\delta dq from dT/dt J/m2/year')

climcaps40 = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_zonal_rates_stats_Sept2002_Aug2021_19yr_desc.mat');
atm_dq_from_dTzdt_per_year = zeros(97,40);
atm_dq_from_dTzdt_per_year(1:97,:)  = mass_a(1:97,:) * 1003.5 .* climcaps40.thestats.ptemprate(:,1:97)';
pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),atm_dq_from_dTzdt_per_year);
set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('\delta dq from dT/dt J/m2/year')
cx = caxis;

figure(2)
%% now suppose 0.005 of this changes per year and then condenses
%% specific heat capacity of air is 1003.5 J/kg/K
atm_heat_from_dWVdt_per_year = mass_w*0.005*2.257e6;
atm_heat_from_dWVdt_per_sec  = mass_w*0.005*2.257e6/(60*60*24*365);
%% actualy we have the WVfrac change
junkclimcapsWV = zeros(size(mass_w)); junkclimcapsWV((1:66)+34-3,:) = climcaps40.thestats.waterrate';
atm_heat_from_dWVdt_per_year = mass_w .* junkclimcapsWV *2.257e6;
atm_heat_from_dWVdt_per_sec  = atm_heat_from_dWVdt_per_year/(60*60*24*365);
pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),atm_heat_from_dWVdt_per_year);
set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('\delta WV condense heat/yr J/m2/year')
caxis(cx);

figure(3);
atm_temp_change_per_sec  = atm_heat_from_dWVdt_per_sec./(mass_a * 1003.5);
atm_temp_change_per_year = atm_heat_from_dWVdt_per_sec./(mass_a * 1003.5) * (60*60*24*365);
pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),atm_temp_change_per_year);
set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(jet); ylim([100 1000]); shading interp; title('\delta WV condense --> dT/dt K/yr')
caxis([-0.15 +0.15])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
coss = ones(97,1)*cos(p.rlat*pi/180);
total_atm_dq_from_dTzdt_per_year   = nansum(atm_dq_from_dTzdt_per_year.*coss,2);
total_atm_heat_from_dWVdt_per_year = nansum(atm_heat_from_dWVdt_per_year.*coss,2);
semilogy(total_atm_dq_from_dTzdt_per_year,p.plays(1:97,20),total_atm_heat_from_dWVdt_per_year,p.plays(1:97,20))
set(gca,'ydir','reverse'); ylim([10 1000]); xlim([-2 +2]*1e5);
plotaxis2;

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load('llsmap5');
figure(1); colormap(llsmap5); caxis([-2 +2]*1e4);
figure(2); colormap(llsmap5); caxis([-2 +2]*1e4);
figure(3); colormap(llsmap5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
plot(p.rlat,climcaps40.thestats.stemprate); plotaxis2;
total_ocean_heat_per_year = climcaps40.thestats.stemprate * 6e6;   %% look at Brian Rose CkimLab notebook, C = 4e8 J/m2 for 100 m layer of ocean, so this is 1.5 m

figure(5);
cossx = ones(97,1)*cos(p.rlat*pi/180);
lat_atm_dq_from_dTzdt_per_year   = nansum(atm_dq_from_dTzdt_per_year.*coss,1);
lat_atm_heat_from_dWVdt_per_year = nansum(atm_heat_from_dWVdt_per_year.*coss,1);
plot(p.rlat,lat_atm_dq_from_dTzdt_per_year/1000,'x-',p.rlat,lat_atm_heat_from_dWVdt_per_year/1000,p.rlat,total_ocean_heat_per_year/1000)
plotaxis2; xlabel('Latitude'); ylabel('kJ/m2/yr heat change');
hl = legend('dT(z)/dt','latent heat dWV(z)/dt','1.5m thick oceans','location','best','fontsize',10);

atmospheric_amplification = ones(97,1) * climcaps40.thestats.stemprate;
atmospheric_amplification = climcaps40.thestats.ptemprate(:,1:97)' ./ atmospheric_amplification;
figure(6); pcolor(ones(97,1)*p.rlat,p.plays(1:97,:),atmospheric_amplification); colorbar;
set(gca,'ydir','reverse'); set(gca,'yscale','log'); colorbar;
colormap(llsmap5); caxis([-2 2]); ylim([10 1000]); shading interp; title('Atmospheric Amplification')
colormap(llsmap5); caxis([-1 1]*3); ylim([10 1000]); shading interp; title('Atmospheric Amplification')

figure(7);
tropics = find(abs(p.rlat) <= 30);
figure(7); plot(atmospheric_amplification(1:97,tropics),p.plays(1:97,tropics),'c',nanmean(atmospheric_amplification(1:97,tropics),2),nanmean(p.plays(1:97,tropics),2),'rx-')
set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); plotaxis2;
