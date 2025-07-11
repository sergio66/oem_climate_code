addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
globalavg = find_zonal_average_rtp(h,p);
if ~isfield(globalavg,'plays')
  globalavg.plays = plevs2plays(globalavg.plevs);
end
[ggLAY,ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2gg(h,globalavg,1:64,1); ggLAY(98:101,:) = NaN; globalavg.wv_gg = ggLAY;
[ggLAY,ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2gg(h,globalavg,1:64,3); ggLAY(98:101,:) = NaN; globalavg.o3_gg = ggLAY;

figure(1); clf; plot(globalavg.rlat,globalavg.stemp)
  xlabel('Latitude'); title('stemp K')

figure(2); clf; pcolor(globalavg.rlat,nanmean(globalavg.plays,2),globalavg.ptemp); colorbar
figure(2); clf; pcolor(globalavg.rlat,globalavg.plays(1:97,:),globalavg.ptemp(1:97,:)); colorbar
figure(2); clf; pcolor(globalavg.rlat,globalavg.plays(1:96,:),globalavg.ptemp(1:96,:)); colorbar
  caxis([200 300]); shading interp; set(gca,'ydir','reverse');
  set(gca,'yscale','log');
  ylim([10 1000])
  xlabel('Latitude'); title('ptemp K')

figure(3); clf; pcolor(globalavg.rlat,nanmean(globalavg.plays,2),globalavg.gas_1); colorbar
figure(3); clf; pcolor(globalavg.rlat,globalavg.plays(1:97,:),globalavg.gas_1(1:97,:)); colorbar
figure(3); clf; pcolor(globalavg.rlat,globalavg.plays(1:96,:),globalavg.gas_1(1:96,:)); colorbar
  caxis([1e15 1e22]); shading interp; set(gca,'ydir','reverse');
  set(gca,'yscale','linear');
  ylim([10 1000])
  xlabel('Latitude'); title('wv molecules/cm2')

figure(4); clf; pcolor(globalavg.rlat,nanmean(globalavg.plays,2),globalavg.gas_3); colorbar
figure(4); clf; pcolor(globalavg.rlat,globalavg.plays(1:97,:),globalavg.gas_3(1:97,:)); colorbar
figure(4); clf; pcolor(globalavg.rlat,globalavg.plays(1:96,:),globalavg.gas_3(1:96,:)); colorbar
  caxis([1e15 1e22]/25000); shading interp; set(gca,'ydir','reverse');
  set(gca,'yscale','log');
  ylim([0.1 1000])
  xlabel('Latitude'); title('o3 molecules/cm2')

figure(5); clf; pcolor(globalavg.rlat,nanmean(globalavg.plays,2),globalavg.wv_gg); colorbar
figure(5); clf; pcolor(globalavg.rlat,globalavg.plays(1:97,:),globalavg.wv_gg(1:97,:)); colorbar
figure(5); clf; pcolor(globalavg.rlat,globalavg.plays(1:96,:),globalavg.wv_gg(1:96,:)); colorbar
  caxis([0 1]/100); shading interp; set(gca,'ydir','reverse');
  set(gca,'yscale','linear');
  ylim([10 1000])
  xlabel('Latitude'); title('wv g/g')

figure(6); clf; pcolor(globalavg.rlat,nanmean(globalavg.plays,2),globalavg.o3_gg); colorbar
figure(6); clf; pcolor(globalavg.rlat,globalavg.plays(1:97,:),globalavg.o3_gg(1:97,:)); colorbar
figure(6); clf; pcolor(globalavg.rlat,globalavg.plays(1:96,:),globalavg.o3_gg(1:96,:)); colorbar
  caxis([0 1]/50000); shading interp; set(gca,'ydir','reverse');
  set(gca,'yscale','log');
  ylim([0.1 1000])
  xlabel('Latitude'); title('o3 g/g')

comment = 'see zonally_averaged_profiles.m';
% save zonally_averaged_profiles.mat h globalavg comment

figure(7); clf; 
junk = globalavg.wv_gg ./ globalavg.gas_1;
pcolor(globalavg.rlat,globalavg.plays(1:96,:),junk(1:96,:)); colorbar
  shading interp; set(gca,'ydir','reverse');
  set(gca,'yscale','linear');
  ylim([100 1000])
  caxis([0 1]*5e-24)
  xlabel('Latitude'); title('wv g/g / molecules/cm2')
