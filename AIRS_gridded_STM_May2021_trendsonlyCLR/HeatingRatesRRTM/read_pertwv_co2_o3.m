addpath /asl/matlib/h4tools;

%{
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm.rtp');

boo = find(p.plevs(:,20) <= 50);
  p.gas_1(boo,:) = p.gas_1(boo,:) .* (1 + 0.001);
boo = find(p.plevs(:,20) > 50);
  p.gas_1(boo,:) = p.gas_1(boo,:) .* (1 + 0.005);

p.gas_2 = p.gas_2 .* (1 + 2.2/385);

boo = find(p.plevs(:,20) >= 1 & p.plevs(:,20) <= 10);
  p.gas_3(boo,:) = p.gas_3(boo,:) .* (1 + 0.01);  %% increase strat O3, from ERA5
boo = find(p.plevs(:,20) >= 10 & p.plevs(:,20) <= 100);
  p.gas_3(boo,:) = p.gas_3(boo,:) .* (1 - 0.01);  %% decrease UT O3, from ERA5

rtpwrite('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/pertWV_CO2_O3.rtp',h,ha,p,pa);
%}


%{
[h,ha,p,pa] = rtpread('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/pertWV_CO2_O3.rtp');
cd /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM
dataFlux = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,p,0);
save /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/PERTWV_CO2_O3/dataFlux_pertWV_CO2_O3.mat dataFlux
mver = ['!/bin/mv /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/DRIVER_CODE_RRTM_Band17/Output_RTPLOOP/00/*.mat /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/PERTWV_CO2_O3/.'];
eval(mver)
%}

cd /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/
for ii = 1 : 40
  fname = ['PERTWV_CO2_O3/rrtm_3p3_plevs_4dp_prof_' num2str(ii) '_emiss0p980.mat'];
  loader = ['a = load(''' fname ''');'];
  eval(loader);
  pertWV_CO2_O3.pres(:,ii) = a.saveinfo.heating_rate_info(1,:,2);
  pertWV_CO2_O3.heatrate(:,ii) = a.saveinfo.heating_rate_info(1,:,6);
end

figure(6); 
pcolor(p.rlat,pertWV_CO2_O3.pres,pertWV_CO2_O3.heatrate); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([0.001 1000]);
colormap jet; xlabel('latitude'); colorbar; ylabel('P (mb)'); shading interp; title('PERTWV\_CO2\_O3 heating rate K/day');

addpath /home/sergio/MATLABCODE/COLORMAP
figure(7); 
pcolor(p.rlat,pertWV_CO2_O3.pres,(pertWV_CO2_O3.heatrate-unpert.heatrate)*365); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([0.001 1000]);
colormap(usa2); caxis([-0.05 +0.05]*365); xlabel('latitude'); colorbar; ylabel('P (mb)'); shading interp; title('PERTWV\_CO2\_O3 heating rate K/yr');

for ii = 1:7;
  figure(ii); ylim([1 1000])
end

figure(3); caxis([-0.05 +0.05]/10*365);
figure(5); caxis([-0.05 +0.05]/10*365);
figure(7); caxis([-0.05 +0.05]/10*365);