addpath /asl/matlib/h4tools;

%{
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm.rtp');
boo = find(p.plevs(:,20) >= 1 & p.plevs(:,20) <= 10);
  p.gas_3(boo,:) = p.gas_3(boo,:) .* (1 + 0.01);  %% increase strat O3, from ERA5
boo = find(p.plevs(:,20) >= 10 & p.plevs(:,20) <= 100);
  p.gas_3(boo,:) = p.gas_3(boo,:) .* (1 - 0.01);  %% decrease UT O3, from ERA5
rtpwrite('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/pertO3.rtp',h,ha,p,pa);
%}


%{
[h,ha,p,pa] = rtpread('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/pertO3.rtp');
cd /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM
dataFlux = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,p,0);
save /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/PERTO3/dataFlux_pertO3.mat dataFlux
mver = ['!/bin/mv /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/DRIVER_CODE_RRTM_Band17/Output_RTPLOOP/00/*.mat /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/PERTO3/.'];
eval(mver)
%}

cd /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/HeatingRatesRRTM/
for ii = 1 : 40
  fname = ['PERTO3/rrtm_3p3_plevs_4dp_prof_' num2str(ii) '_emiss0p980.mat'];
  loader = ['a = load(''' fname ''');'];
  eval(loader);
  pertO3.pres(:,ii) = a.saveinfo.heating_rate_info(1,:,2);
  pertO3.heatrate(:,ii) = a.saveinfo.heating_rate_info(1,:,6);
end

figure(4); 
pcolor(p.rlat,pertO3.pres,pertO3.heatrate); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([0.001 1000]);
colormap jet; xlabel('latitude'); colorbar; ylabel('P (mb)'); shading interp; title('PERTO3 heating rate K/day');

addpath /home/sergio/MATLABCODE/COLORMAP
figure(5); 
pcolor(p.rlat,pertO3.pres,(pertO3.heatrate-unpert.heatrate)*365); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([0.001 1000]);
colormap(usa2); caxis([-0.05 +0.05]*365); xlabel('latitude'); colorbar; ylabel('P (mb)'); shading interp; title('PERTO3 heating rate K/yr');

for ii = 1:5;
  figure(ii); ylim([1 1000])
end
