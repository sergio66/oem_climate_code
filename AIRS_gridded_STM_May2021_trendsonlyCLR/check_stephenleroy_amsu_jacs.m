iIRjac = input('do you want to make SARTA AIRS jacs (-1 default/+1 : ')
iIRjac = -1;

if length(iIRjac) == 0
  iIRjac = -1;
end
if iIRjac > 0
  sartaIR = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 '];
  sartaIR = [sartaIR ' fin=summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp fout=junk.jac.rtp'];
  sartaIR = [sartaIR ' listp=' num2str(driver.iibin) ' listj=100,1 jacunit=2'];
  eval(sartaIR)
  addpath /home/sergio/KCARTA/MATLAB
  figure(9); clf
  figure(10); clf
  [sartaW,sartaJACTz,~,iaNumLay] = readsarta_jacV2('junk.jac.rtp_jacTZ',100); sartaJACTz = sartaJACTz';
  [sartaW,sartaJACG1,~,iaNumLay] = readsarta_jacV2('junk.jac.rtp_jacG1',001); sartaJACG1 = sartaJACG1';
  eval(['!/bin/rm junk.jac.rtp_jacTZ junk.jac.rtp_jacG1 junk.jac.rtp']);
  sartaJACTz = squeeze(sartaJACTz); sartaJACSKT = sartaJACTz(:,101); sartaJACTz = sartaJACTz(:,1:100);
  sartaJACG1 = squeeze(sartaJACG1);
  figure(09); pcolor(sartaW,1:97,sartaJACTz(:,1:97)'); 
    shading interp; set(gca,'ydir','reverse'); colorbar; colormap(usa2); xlim([640 1640]); title('IR Tz jac')
  figure(10); pcolor(sartaW,1:97,sartaJACG1(:,1:97)'); 
    shading interp; set(gca,'ydir','reverse'); colorbar; colormap(usa2); xlim([640 1640]); title('IR WV jac')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
junk = read_netcdf_lls('../OSS_AMSU_jacs/oss_bt_jac_avg.cloudy.2010.nc');
bt_amsu = squeeze(nanmean(squeeze(nanmean(junk.brightness_temperature,4)),3));
figure(1); clf; plot(junk.lats,bt_amsu(1,:))
figure(1); clf; plot(junk.lats,bt_amsu(1:2,:),'+-',junk.lats,bt_amsu(3:11,:))
  legend('Ch4','Ch5','location','best'); ylabel('Men BT'); xlabel('Latitude');

sleroy.T    = squeeze(nanmean(junk.temperature_jacobian(:,:,driver.iLat,driver.iLon,:),5));
sleroy.WV   = squeeze(nanmean(junk.humidity_jacobian(:,:,driver.iLat,driver.iLon,:),5));
sleroy.P    = squeeze(nanmean(junk.levels(:,driver.iLat,driver.iLon,:),4));
sleroy.SKT  = squeeze(nanmean(junk.tskin_jacobian(:,driver.iLat,driver.iLon,:),4));
sleroy.emis = squeeze(nanmean(junk.emissivity_jacobian(:,driver.iLat,driver.iLon,:),4));

figure(1); 
  pcolor([jac.f(4:end);57.29],sleroy.P,sleroy.T'); colorbar;
  colormap jet; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]);
  xlabel('AMSU channel freq GHz'); ylabel('P [mb]');
  title('T jac tile S. Leroy')
  ax1 = axis;

figure(2); 
  pcolor([jac.f(4:end);57.29],sleroy.P,sleroy.WV'); colorbar;
  colormap jet; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]);
  xlabel('AMSU channel freq GHz'); ylabel('P [mb]');
  title('WV jac tile S. Leroy')
  ax2 = axis;
  caxis([-1/3 +1]*30)

playsx = load('/home/sergio/MATLABCODE/airslevels.dat');
playsx = flipud(plevs2plays(playsx));
playsx = playsx(5:101);

%% /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/driver_amsu_jacs.m
log10 = log(10);  %% did I use this, nah, see driver_amsu_jacs.m and clust_driver_AMSU_jacs.m
log10 = 1;

figure(3); 
  pcolor(jac.f,playsx,m_ts_jac0(:,driver.jacobian.temp_i)'); colorbar;
  colormap jet; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]);
  xlabel('AMSU channel freq GHz'); ylabel('P [mb]');
  title('T jac tile SARTA P.Rosenkranz')
  axis(ax1);

figure(4); 
  hmm = m_ts_jac0(:,driver.jacobian.water_i)';
  pcolor(jac.f,playsx,m_ts_jac0(:,driver.jacobian.water_i)'); colorbar;
  colormap jet; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]);
  xlabel('AMSU channel freq GHz'); ylabel('P [mb]');
  title('WV jac tile SARTA P.Rosenkranz')
  caxis([0 1]/50)
  axis(ax2);

aha = load('zonally_averaged_profiles.mat');
ratio = aha.globalavg.wv_gg ./ aha.globalavg.gas_1;
ratio = log(ratio);

ratio = flipud(flipud(1./aha.globalavg.wv_gg));
ratio = 1./aha.globalavg.wv_gg;

ratio = ratio(:,driver.iLat);
ratio_save = ratio;

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8); loglog(ratio_save(5:101),playsx); set(gca,'ydir','reverse'); ylabel('P[mb]'); xlabel('1/SH'); ylim([50 1000])
figure(8); semilogy(hmm(:,4),playsx); set(gca,'ydir','reverse');ylabel('P[mb]');  xlabel('WV jac'); ylim([50 1000])

x1 = ratio_save(5:101);
y1 = playsx;
x2 = hmm(:,4);
y2 = playsx;

figure(8); clf
t = tiledlayout(1,1);
ax1 = axes(t);
loglog(ax1,x1,y1,'-r');
ax1.XColor = 'r';
ax1.YColor = 'r';
set(gca,'ydir','reverse'); ylim([50 1000]); xlabel('1/SH');
%
ax2 = axes(t);
semilogy(ax2,x2,y2,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
set(gca,'ydir','reverse'); ylim([50 1000]); xlabel('WV jac');

%%%%%%%%%%%%%%%%%%%%%%%%%

ratio = ratio(5:101) * ones(1,length(jac.f));;
ratio = ratio';

figure(5); 
  pcolor(jac.f,playsx,log10*(ratio.*m_ts_jac0(:,driver.jacobian.water_i))'); colorbar;
  colormap jet; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]);
  xlabel('AMSU channel freq GHz'); ylabel('P [mb]');
  title('WV jac tile SARTA P.Rosenkranz')
  axis(ax2);
  caxis([-1/3 +1]*3)

figure(6); clf
  plot([jac.f(4:end);57.29],sleroy.SKT,'b',jac.f,m_ts_jac0(:,6),'c',...
       [jac.f(4:end);57.29],nansum(sleroy.T,2),'r',jac.f,nansum(m_ts_jac0(:,driver.jacobian.temp_i),2),'m','linewidth',2)
xlim([ax2(1) ax2(2)])
legend('SL SKT','SSM SKT','SL sum(T)','SSM sum(T)','location','best');

figure(7);
  plot([jac.f(4:end);57.29],nansum(sleroy.WV,2),'g',jac.f,nansum(ratio.*m_ts_jac0(:,driver.jacobian.water_i),2)*log10,'k','linewidth',2)
xlim([ax2(1) ax2(2)])
legend('SL sum(WV)','SSM sum(WV)/q','location','best');
