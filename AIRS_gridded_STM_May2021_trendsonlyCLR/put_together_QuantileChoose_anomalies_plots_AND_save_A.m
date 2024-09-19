figure(1); clf; pcolor(reshape(bt1231_A,72,64)'); title('BT 1231'); shading interp; colorbar; xlabel('LonBin 01-72'); ylabel('LatBin 01-64'); colormap jet; pause(0.1);
figure(2); clf; pcolor(reshape(btChID_A,72,64)'); title('BT ChID'); shading interp; colorbar; xlabel('LonBin 01-72'); ylabel('LatBin 01-64'); colormap jet; pause(0.1);
datime = 2002.75 + (1:iNumAnomTimeSteps)*16/365;
  figure(3); clf; junkjunk = squeeze(btanomA(:,i1231,:));  junkjunk = reshape(junkjunk,72,64,iNumAnomTimeSteps); junkjunk = squeeze(nanmean(junkjunk,1)); 
    pcolor(datime,1:64,junkjunk); pcolor(datime,1:64,smoothn(junkjunk,1)); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*2); title('BT 1231 anom');
  figure(4); clf; junkjunk = squeeze(btanomA(:,iXCh,:)); junkjunk = reshape(junkjunk,72,64,iNumAnomTimeSteps); junkjunk = squeeze(nanmean(junkjunk,1)); 
    pcolor(datime,1:64,junkjunk); pcolor(datime,1:64,smoothn(junkjunk,1)); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*2); title('BT ChID anom');
pause(0.1) 

figure(5); clf; pcolor(reshape(bt1231_D,72,64)'); title('desc BT 1231'); colormap jet; shading interp; colorbar; ; caxis([200 320])
figure(6); clf; pcolor(reshape(btChID_D,72,64)'); title('desc BT ChID'); colormap jet; shading interp; colorbar; ;
figure(7); clf; pcolor(reshape(bt1231_D-btChID_D,72,64)'); title('desc BT1231-BTXXXX'); colormap(usa2); shading interp; colorbar; ; caxis([-1 +1]*5)
figure(8); clf; pcolor(reshape(bt1231_A,72,64)'); title('asc BT 1231'); colormap jet; shading interp; colorbar; ;  caxis([200 320])
figure(9); clf; pcolor(reshape(btChID_A,72,64)'); title('asc BT ChID'); colormap jet; shading interp; colorbar; ;
figure(10); clf; pcolor(reshape(bt1231_A-btChID_A,72,64)'); title('assc BT1231-BTXXXX'); colormap(usa2); shading interp; colorbar; ; caxis([-1 +1]*5)
pause(0.1);

a1 = load(fname,'yy_desc','mm_desc','dd_desc','hh_desc','rtime_desc');
yyD = a1.yy_desc;
mmD = a1.mm_desc;
ddD = a1.dd_desc;
hhD = a1.hh_desc;
rtimeD = a1.rtime_desc;

a1 = load(fname,'yy_asc','mm_asc','dd_asc','hh_asc','rtime_asc');
yyA = a1.yy_asc;
mmA = a1.mm_asc;
ddA = a1.dd_asc;
hhA = a1.hh_asc;
rtimeA = a1.rtime_asc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/COLORMAP

do_XX_YY_from_X_Y
[salti, landfrac] = usgs_deg10_dem(Y,X);
lf = landfrac; lf = lf(:);
[mmjunk,nnjunk] = size(XX);
if mmjunk == 1
  XX = XX';
  YY = YY';
end

iOffSet = 4;

figure(iOffSet+1); clf; scatter_coast(XX(:),YY(:),100,bt1231_D(:)); title('desc BT 1231'); colormap jet; shading interp; colorbar; ; caxis([200 320])
figure(iOffSet+2); clf; scatter_coast(XX(:),YY(:),100,btChID_D(:)); title('desc BT ChID'); colormap jet; shading interp; colorbar; ;
figure(iOffSet+3); clf; scatter_coast(XX(:),YY(:),100,bt1231_D(:)-btChID_D(:)); title('desc BT1231-BTXXXX'); colormap(usa2); shading interp; colorbar; ; caxis([-1 +1]*5)
figure(iOffSet+4); clf; scatter_coast(XX(:),YY(:),100,bt1231_A(:)); title('asc BT 1231'); colormap jet; shading interp; colorbar; ;  caxis([200 320])
figure(iOffSet+5); clf; scatter_coast(XX(:),YY(:),100,btChID_A(:)); title('asc BT ChID'); colormap jet; shading interp; colorbar; ;
figure(iOffSet+6); clf; scatter_coast(XX(:),YY(:),100,bt1231_A(:)-btChID_A(:)); title('assc BT1231-BTXXXX'); colormap(usa2); shading interp; colorbar; ; caxis([-1 +1]*5)
pause(0.1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coslat = cos(YY*pi/180)*ones(1,iNumAnomTimeSteps);

if iHuge < 0
  figure(iOffSet+7); clf
   daysSince2002A = change2days(yyA,mmA,ddA,2002);
   P = polyfit(daysSince2002A,nanmean(btanomA,1),1);         fprintf(1,'quick BTtrend = %8.6f K/year \n',P(1)*365)
   P = polyfit(daysSince2002A,nanmean(btanomA.*coslat,1),1); fprintf(1,'wgted BTtrend = %8.6f K/year \n',P(1)*365)
   
   pcolor(2002+daysSince2002A/365,rlat,squeeze(nanmean(reshape(btanomA,72,64,iNumAnomTimeSteps),1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1)
     xlabel('Time'); ylabel('Latitude'); title(['BT' num2str(h.vchan(iChIA),'%6.2f') ' cm-1 anomaly']); 
   hold on; plot(2002+daysSince2002A/365,100*nanmean(btanomA,1),'k','linewidth',2); hold off; plotaxis2;
   
  figure(iOffSet+8); clf
   oni = load('ONI_sep2023.txt');
   [aaa,bbb] = size(oni);
   oniS = oni(1,1); oniE = oni(aaa,1);
   oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
   onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
   plot(onidd,10*oni,'k',2002+daysSince2002A/365,smooth(100*nanmean(btanomA,1),10),'b','linewidth',2); plotaxis2; xlim([2002 2023]); ylim([-1 +1]*50)
   
  P = polyfit(daysSince2002A,nanmean(btanomA,1),1); YP = polyval(P,daysSince2002A);
  PX = polyfit(daysSince2002A,nanmean(btanomA.*coslat,1),1); YPX = polyval(PX,daysSince2002A);
  
  plot(onidd,10*oni,'k',...
       2002+daysSince2002A/365,smooth(100*(nanmean(btanomA,1)),10),'b',...
       2002+daysSince2002A/365,smooth(100*(nanmean(btanomA,1)-YP),10),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002A/365,smooth(100*(nanmean(btanomA,1)-YP),24),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002A/365,smooth(100*(nanmean(btanomA.*coslat,1)-YPX),24),'r','linewidth',2); xlim([2002 2023]);
  plotaxis2;
  hl = legend('ONI',['BT' num2str(h.vchan(iChIA),'%6.2f') ' cm-1 anomaly'],'location','best');
  hl = legend('ONI',['BT' num2str(h.vchan(iChIA),'%6.2f')],'location','best');

  [U,S,V] = svd(btanomA,'econ');

elseif iHuge > 0

  figure(iOffSet+7);
  daysSince2002A = change2days(yyA,mmA,ddA,2002);
  daysSince2002D = change2days(yyD,mmD,ddD,2002);
  moo = squeeze(nanmean(btanomA,1));
  P = polyfit(daysSince2002A,moo(i1231,:),1);             fprintf(1,'quick BTtrend = %8.6f K/year \n',P(1)*365)
  moox = squeeze(btanomA(:,i1231,:));
  P = polyfit(daysSince2002A,nanmean(moox.*coslat,1),1); fprintf(1,'wgted BTtrend = %8.6f K/year \n',P(1)*365)
  
  pcolor(2002+daysSince2002A/365,rlat,squeeze(nanmean(reshape(moox,72,64,iNumAnomTimeSteps),1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1)
    xlabel('Time'); ylabel('Latitude'); title(['BT' num2str(h.vchan(iChIA),'%6.2f') ' cm-1 anomaly']); 
  hold on; plot(2002+daysSince2002A/365,100*nanmean(moox,1),'k','linewidth',2); hold off; plotaxis2;
  
  figure(iOffSet+8);
  oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
  plot(onidd,10*oni,'k',2002+daysSince2002A/365,smooth(100*nanmean(moox,1),10),'b','linewidth',2); xlim([2002 2023]);
  
  P = polyfit(daysSince2002A,nanmean(moox,1),1); YP = polyval(P,daysSince2002A);
  PX = polyfit(daysSince2002A,nanmean(moox.*coslat,1),1); YPX = polyval(PX,daysSince2002A);
  
  plot(onidd,10*oni,'k',...
       2002+daysSince2002A/365,smooth(100*(nanmean(moox,1)),10),'b',...
       2002+daysSince2002A/365,smooth(100*(nanmean(moox,1)-YP),10),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002A/365,smooth(100*(nanmean(moox,1)-YP),24),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002A/365,smooth(100*(nanmean(moox.*coslat,1)-YPX),24),'r','linewidth',2); xlim([2002 2023]);
  plotaxis2;
  hl = legend('ONI',['BT' num2str(h.vchan(iChIA),'%6.2f') ' cm-1 anomaly'],'location','best');
  hl = legend('ONI',['BT' num2str(h.vchan(iChIA),'%6.2f')],'location','best');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(iOffSet+9); clf

%% raw data in eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin38/LonBin51/iQAX_3_summarystats_LatBin38_LonBin51_timesetps_001_457_V1.mat
if iHuge < 0
  if exist('btanonA')
    plot(2002 + daysSince2002A/365,btanomA(2787,:),'b',2002 + daysSince2002A/365,btanomA(2787,:),'r','linewidth',2); %% LatBin = 39, LonBin = 51 >> (39-1)*72 + 51 = 2787
  else
    plot(2002 + daysSince2002A/365,btanomA(2787,:)*0,'b',2002 + daysSince2002A/365,btanomA(2787,:),'r','linewidth',2); %% LatBin = 39, LonBin = 51 >> (39-1)*72 + 51 = 2787
  end
  plotaxis2;
  hl = legend('asc','desc','location','best'); title(['BTXXXX anom, tile 2787 over India/Arabian Sea Q' num2str(iQuant)]); xlim([2002.75 2022.75])
end

%%%%%%%%%%%%%%%%%%%%%%%%%

%% raw data in eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin38/LonBin52/iQAX_3_summarystats_LatBin38_LonBin52_timesetps_001_457_V1.mat
%{
%
% moo = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin38/LonBin52/iQAX_3_summarystats_LatBin38_LonBin52_timesetps_001_457_V1.mat');
% mootime = moo.year_desc + (moo.month_desc-1)/12 + (moo.day_desc-1)/12/30;
% 
% plot(mootime,moo.meanBT1231_asc+10,'b',mootime,moo.meanBT1231_desc-10,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile 2788')
% B = Math_tsfit_lin_robust((mootime-2002)*365,moo.meanBT1231_asc,4);  fprintf(1,'BT1231 trend for asc node = %8.6f K/yr \n',B(2))
% B = Math_tsfit_lin_robust((mootime-2002)*365,moo.meanBT1231_desc,4); fprintf(1,'BT1231 trend for desc node = %8.6f K/yr \n',B(2))
% 
% moo.rad1231_asc  = squeeze(moo.rad_quantile_asc(:,i1231,5));
% moo.rad1231_desc = squeeze(moo.rad_quantile_desc(:,i1231,5));
% plot(mootime,moo.rad1231_asc+0,'b',mootime,moo.rad1231_desc-0,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile 2788, Q05')
% plot(mootime,moo.rad1231_asc+5,'b',mootime,moo.rad1231_desc-5,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile 2788, Q05')
% B = Math_tsfit_lin_robust((mootime-2002)*365,moo.rad1231_asc,4);  fprintf(1,'rad1231 trend Q05 for asc node = %8.6f K/yr \n',B(2))
% B = Math_tsfit_lin_robust((mootime-2002)*365,moo.rad1231_desc,4); fprintf(1,'rad1231 trend Q05 for desc node = %8.6f K/yr \n',B(2))
% 
% moo.bt1231_asc  = rad2bt(1231,squeeze(moo.rad_quantile_asc(:,i1231,5)));
% moo.bt1231_desc = rad2bt(1231,squeeze(moo.rad_quantile_desc(:,i1231,5)));
% plot(mootime,moo.bt1231_asc+0,'b',mootime,moo.bt1231_desc-0,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile 2788, Q05')
% %plot(mootime,moo.bt1231_asc+5,'b',mootime,moo.bt1231_desc-5,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile 2788, Q05')
% B = Math_tsfit_lin_robust((mootime-2002)*365,moo.bt1231_asc,4);  fprintf(1,'bt1231 trend Q05 for asc node = %8.6f K/yr \n',B(2))
% B = Math_tsfit_lin_robust((mootime-2002)*365,moo.bt1231_desc,4); fprintf(1,'bt1231 trend Q05 for desc node = %8.6f K/yr \n',B(2))
% 
% sktjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/surface_jac_new.mat');
% tjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/temp_jac_new.mat');
% g1jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac_new.mat');
% g101jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g101_jac_new.mat');
% g102jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac_new.mat');
% plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2))
% plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2)*0.1); xlim([640 1640]); plotaxis2;
% 
% sktjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/surface_jac.mat');
% tjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/temp_jac.mat');
% g1jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g1_jac.mat');
% g101jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g101_jac.mat');
% g102jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g1_jac.mat');
% plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2))
% plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2)*0.1); xlim([640 1640]); plotaxis2;
% 
%}

if iHuge < 0
  figure(iOffSet+9); clf
  if exist('btanomA')
    plot(2002 + daysSince2002A/365,btanomA(2788,:),'b',2002 + daysSince2002A/365,btanomA(2788,:),'r','linewidth',2); %% LatBin = 39, LonBin = 52 >> (39-1)*72 + 52 = 2788
  else
    plot(2002 + daysSince2002A/365,btanomA(2788,:)*0,'b',2002 + daysSince2002A/365,btanomA(2788,:),'r','linewidth',2); %% LatBin = 39, LonBin = 52 >> (39-1)*72 + 52 = 2788
  end
  plotaxis2;
  hl = legend('asc','desc','location','best'); title(['BTXXXX anom, tile 2788 over Central India Q' num2str(iQuant)]); xlim([2002.75 2022.75])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSave = input('save (-1/+1) : ');
if iSave > 0
  comment = 'look at driver_put_together_QuantileChoose_anomalies.m';
  if iHuge < 0
    fout = ['anomaly_chID_' num2str(iChID,'%04d') '_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat'];  %% orig
    fout = ['anomaly_iQAX_' num2str(iQAX) '_chID_' num2str(iChID,'%04d') '_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat'];
    saver = ['save ' fout ' bt1231_D btChID_D btanomD yyD mmD ddD hhD rtimeD fluxanomD    btanomA yyA mmA ddA hhA rtimeA bt1231_A btChIDn_A  fluxanomA  comment ind f_ind'];
  else
    fout = ['/asl/s1/sergio/JUNK/anomaly_ALL_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '_A.mat'];  %% orig
    fout = ['/asl/s1/sergio/JUNK/anomaly_iQAX_' num2str(iQAX) '_ALL_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '_A.mat'];
    saver = ['save -v7.3 ' fout ' btanomA yyA mmA ddA hhA rtimeA comment bt1231_A fluxanomA ind f_ind'];
  end
  eval(saver)
end

