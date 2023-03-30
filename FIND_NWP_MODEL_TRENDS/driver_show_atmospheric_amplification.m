pert = make_rtp_plays(pert);
pavg = nanmean(pert.plays,2); boo = find(pert.spres == max(pert.spres)); pavg = pert.plays(:,boo);

disp('Enter (-7) polar    L/O')
disp('      (-6) midlat   L/O')
disp('      (-5) tropical L/O')
disp('      (-4) polar land    (+4) polar ocean')
disp('      (-3) midlat land   (+3) midlat ocean')
disp('      (-2) tropical land (+2) tropical ocean')
disp('      (-1) land          (+1) ocean');
disp('      [0,default] ALL trends : ');
iAorOorL = input('Enter region : ');
if length(iAorOorL) == 0 
  iAorOorL = 0;
end

%% this is bascially copied from AIRS_gridded_STM_May2021_trendsonlyCLR/plot_driver_gather_gridded_retrieval_results. but with TWO modifications
clear maskLF
maskLF = zeros(1,4608);
maskLF = nan(1,4608);    %% MODIFICATION 1
if size(landfrac) ~= size(Xlon)
  Xlon = Xlon';
  Ylat = Ylat';
end
if iAorOorL == -7
  maskLF(abs(Ylat) > 60) = 1;
elseif iAorOorL == -6
  maskLF(abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == -5
  maskLF(abs(Ylat) < 30) = 1;
elseif iAorOorL == 0
  maskLF = ones(1,4608);
elseif iAorOorL == -1
  maskLF(landfrac == 1) = 1;
elseif iAorOorL == +1
  maskLF(landfrac == 0) = 1;
elseif iAorOorL == -2
  maskLF(landfrac == 1 & abs(Ylat) <= 30) = 1;
elseif iAorOorL == +2
  maskLF(landfrac == 0 & abs(Ylat) <= 30) = 1;
elseif iAorOorL == -3
  maskLF(landfrac == 1 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == +3
  maskLF(landfrac == 0 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == -4
  maskLF(landfrac == 1 & abs(Ylat) > 60) = 1;
elseif iAorOorL == +4
  maskLF(landfrac == 0 & abs(Ylat) > 60) = 1;
end
maskLFmatr = reshape(maskLF,72,64)';
mask = find(maskLF == 1);
maskLF100 = ones(100,1)*maskLF;  %% MODIFICATION 2

clear junk;
  junk.pavg = pavg;
  junk.stemprate = reshape(airsL3.stemprate,72,64);     junk.stemprate = reshape(junk.stemprate,1,4608).*maskLF;
  junk.RHrate    = permute(airsL3.RHrate,[3 1 2]);      junk.RHrate    = reshape(junk.RHrate,100,4608).*maskLF100;
  junk.ptemprate = permute(airsL3.ptemprate,[3 1 2]);   junk.ptemprate = reshape(junk.ptemprate,100,4608).*maskLF100;
  junk.waterrate = permute(airsL3.waterrate,[3 1 2]);   junk.waterrate = reshape(junk.waterrate,100,4608).*maskLF100;
     airsL3amp = atmospheric_amplification(junk);
  junkAIRSL3 = junk;

clear junk;
  junk.pavg = pavg;
  junk.stemprate = reshape(climcapsL3.stemprate,72,64);     junk.stemprate = reshape(junk.stemprate,1,4608).*maskLF;
  junk.RHrate    = permute(climcapsL3.RHrate,[3 1 2]);      junk.RHrate    = reshape(junk.RHrate,100,4608).*maskLF100;
  junk.ptemprate = permute(climcapsL3.ptemprate,[3 1 2]);   junk.ptemprate = reshape(junk.ptemprate,100,4608).*maskLF100;
  junk.waterrate = permute(climcapsL3.waterrate,[3 1 2]);   junk.waterrate = reshape(junk.waterrate,100,4608).*maskLF100;
     climcapsL3amp = atmospheric_amplification(junk);
  junkCLIMCAPSL3 = junk;

clear junk
  junk = era5;
  junk.pavg = pavg;
     junk.stemprate = junk.stemprate .* maskLF;
     junk.ptemprate = junk.ptemprate .* maskLF100;
     junk.waterrate = junk.waterrate .* maskLF100;
     junk.RHrate = junk.RHrate .* maskLF100;
     era5amp = atmospheric_amplification(junk);
  junkERA5 = junk;

clear junk
  junk = merra2;
  junk.pavg = pavg;
     junk.stemprate = junk.stemprate .* maskLF;
     junk.ptemprate = junk.ptemprate .* maskLF100;
     junk.waterrate = junk.waterrate .* maskLF100;
     junk.RHrate = junk.RHrate .* maskLF100;
     merra2amp = atmospheric_amplification(junk);
  junkMERRA2 = junk;

iFig = iFigiAmp;
iFig = iFig + 1;
figure(iFig); clf; plot(era5amp.T_amp,pavg(1:97),'b',merra2amp.T_amp,pavg(1:97),'c',airsL3amp.T_amp,pavg(1:97),'r',climcapsL3amp.T_amp,pavg(1:97),'m','linewidth',2)
  set(gca,'ydir','reverse'); plotaxis2;
  title('T amplification dT/dST'); ylim([100 1000]); xlim([-1 +1]); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','location','best','fontsize',10);

iFig = iFig + 1;
figure(iFig); clf; plot(era5amp.RH_amp,pavg(1:97),'b',merra2amp.RH_amp,pavg(1:97),'c',airsL3amp.RH_amp,pavg(1:97),'r',climcapsL3amp.RH_amp,pavg(1:97),'m','linewidth',2)
  set(gca,'ydir','reverse'); plotaxis2;
  title('RH amplification dRH/dST'); ylim([100 1000]); xlim([-1 +1]*5); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','location','best','fontsize',10);

iFig = iFig + 1;
figure(iFig); clf; plot(era5amp.WVfrac_percent_amp,pavg(1:97),'b',merra2amp.WVfrac_percent_amp,pavg(1:97),'c',...
                     airsL3amp.WVfrac_percent_amp,pavg(1:97),'r',climcapsL3amp.WVfrac_percent_amp,pavg(1:97),'m','linewidth',2)
  set(gca,'ydir','reverse'); plotaxis2;
  title('WVfrac amplification dpercentWV/dST'); ylim([100 1000]); xlim([-1 +1]*50/2); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','location','best','fontsize',10);

if iUMBC > 0
  clear junk
  junk = umbc;
  junk.pavg = pavg;
     junk.stemprate = junk.stemprate .* maskLF;
     junk.ptemprate = junk.ptemprate .* maskLF100;
     junk.waterrate = junk.waterrate .* maskLF100;
     junk.RHrate = junk.RHrate .* maskLF100;
     umbcamp = atmospheric_amplification(junk);
  junkUMBC = junk;

  iFig = iFigiAmp;
  iFig = iFig + 1;
  figure(iFig); clf; 
    subplot(121);
      plot(era5amp.T_amp,pavg(1:97),'b',merra2amp.T_amp,pavg(1:97),'c',airsL3amp.T_amp,pavg(1:97),'r',climcapsL3amp.T_amp,pavg(1:97),'m',... 
           umbcamp.T_amp,pavg(1:97),'k','linewidth',2)
      set(gca,'ydir','reverse'); plotaxis2;
      title('T amplification dT/dST'); ylim([100 1000]); xlim([-1 +1]); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);
    subplot(122);
      plot(era5amp.T_amp_sig,pavg(1:97),'b',merra2amp.T_amp_sig,pavg(1:97),'c',airsL3amp.T_amp_sig,pavg(1:97),'r',climcapsL3amp.T_amp_sig,pavg(1:97),'m',... 
           umbcamp.T_amp_sig,pavg(1:97),'k','linewidth',2)
      set(gca,'ydir','reverse'); plotaxis2;
      title('T amplification dT/dST'); ylim([100 1000]); xlim([0 +1]); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);
  
  iFig = iFig + 1;
  figure(iFig); clf; 
    subplot(121);
      plot(era5amp.RH_amp,pavg(1:97),'b',merra2amp.RH_amp,pavg(1:97),'c',airsL3amp.RH_amp,pavg(1:97),'r',climcapsL3amp.RH_amp,pavg(1:97),'m',...
           umbcamp.RH_amp,pavg(1:97),'k','linewidth',2)
      set(gca,'ydir','reverse'); plotaxis2;
      title('RH amplification dRH/dST'); ylim([100 1000]); xlim([-1 +1]*5); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);
    subplot(122);
      plot(era5amp.RH_amp_sig,pavg(1:97),'b',merra2amp.RH_amp_sig,pavg(1:97),'c',airsL3amp.RH_amp_sig,pavg(1:97),'r',climcapsL3amp.RH_amp_sig,pavg(1:97),'m',...
           umbcamp.RH_amp_sig,pavg(1:97),'k','linewidth',2)
      set(gca,'ydir','reverse'); plotaxis2;
      title('RH amplification dRH/dST'); ylim([100 1000]); xlim([0 +1]*5); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);
  
  iFig = iFig + 1;
  figure(iFig); clf; 
    subplot(121);
      plot(era5amp.WVfrac_percent_amp,pavg(1:97),'b',merra2amp.WVfrac_percent_amp,pavg(1:97),'c',...
           airsL3amp.WVfrac_percent_amp,pavg(1:97),'r',climcapsL3amp.WVfrac_percent_amp,pavg(1:97),'m',...
           umbcamp.WVfrac_percent_amp,pavg(1:97),'k','linewidth',2)
      set(gca,'ydir','reverse'); plotaxis2;
      title('WVfrac amplification dpercentWV/dST'); ylim([100 1000]); xlim([-1 +1]*50/2); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);
    subplot(122);
      plot(era5amp.WVfrac_percent_amp_sig,pavg(1:97),'b',merra2amp.WVfrac_percent_amp_sig,pavg(1:97),'c',...
           airsL3amp.WVfrac_percent_amp_sig,pavg(1:97),'r',climcapsL3amp.WVfrac_percent_amp_sig,pavg(1:97),'m',...
           umbcamp.WVfrac_percent_amp_sig,pavg(1:97),'k','linewidth',2)
      set(gca,'ydir','reverse'); plotaxis2;
      title('WVfrac amplification dpercentWV/dST'); ylim([100 1000]); xlim([0 +1]*50/2); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);

  %%%%%%%%%%
  iFig = iFigiAmp;

  iFig = iFig + 1;
  figure(iFig); clf; 
    plot(era5amp.T_amp,pavg(1:97),'b',merra2amp.T_amp,pavg(1:97),'c',airsL3amp.T_amp,pavg(1:97),'r',climcapsL3amp.T_amp,pavg(1:97),'m',... 
           umbcamp.T_amp,pavg(1:97),'k','linewidth',2)
    hold on; 
    shadedErrorBarY(era5amp.T_amp,pavg(1:97),era5amp.T_amp_sig,'b',0.25);
    shadedErrorBarY(merra2amp.T_amp,pavg(1:97),merra2amp.T_amp_sig,'c',0.1);
    shadedErrorBarY(airsL3amp.T_amp,pavg(1:97),airsL3amp.T_amp_sig,'r',0.1);
    shadedErrorBarY(climcapsL3amp.T_amp,pavg(1:97),climcapsL3amp.T_amp_sig,'m',0.1);
    shadedErrorBarY(umbcamp.T_amp,pavg(1:97),umbcamp.T_amp_sig,'k',0.25);
    hold off
    set(gca,'ydir','reverse'); plotaxis2;
    title('T amplification dT/dST'); ylim([100 1000]); xlim([-1 +1]*2); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);

  iFig = iFig + 1;
  figure(iFig); clf; 
    plot(era5amp.RH_amp,pavg(1:97),'b',merra2amp.RH_amp,pavg(1:97),'c',airsL3amp.RH_amp,pavg(1:97),'r',climcapsL3amp.RH_amp,pavg(1:97),'m',... 
           umbcamp.RH_amp,pavg(1:97),'k','linewidth',2)
    hold on; 
    shadedErrorBarY(era5amp.RH_amp,pavg(1:97),era5amp.RH_amp_sig,'b',0.25);
    shadedErrorBarY(merra2amp.RH_amp,pavg(1:97),merra2amp.RH_amp_sig,'c',0.1);
    shadedErrorBarY(airsL3amp.RH_amp,pavg(1:97),airsL3amp.RH_amp_sig,'r',0.1);
    shadedErrorBarY(climcapsL3amp.RH_amp,pavg(1:97),climcapsL3amp.RH_amp_sig,'m',0.1);
    shadedErrorBarY(umbcamp.RH_amp,pavg(1:97),umbcamp.RH_amp_sig,'k',0.25);
    hold off
    set(gca,'ydir','reverse'); plotaxis2;
    title('RH amplification dRH/dST'); ylim([100 1000]); xlim([-1 +1]*5); hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);

  iFig = iFig + 1;
  figure(iFig); clf; 
    plot(era5amp.WVfrac_percent_amp,pavg(1:97),'b',merra2amp.WVfrac_percent_amp,pavg(1:97),'c',...
           airsL3amp.WVfrac_percent_amp,pavg(1:97),'r',climcapsL3amp.WVfrac_percent_amp,pavg(1:97),'m',... 
           umbcamp.WVfrac_percent_amp,pavg(1:97),'k','linewidth',2)
    hold on; 
    shadedErrorBarY(era5amp.WVfrac_percent_amp,pavg(1:97),era5amp.WVfrac_percent_amp_sig,'b',0.25);
    shadedErrorBarY(merra2amp.WVfrac_percent_amp,pavg(1:97),merra2amp.WVfrac_percent_amp_sig,'c',0.1);
    shadedErrorBarY(airsL3amp.WVfrac_percent_amp,pavg(1:97),airsL3amp.WVfrac_percent_amp_sig,'r',0.1);
    shadedErrorBarY(climcapsL3amp.WVfrac_percent_amp,pavg(1:97),climcapsL3amp.WVfrac_percent_amp_sig,'m',0.1);
    shadedErrorBarY(umbcamp.WVfrac_percent_amp,pavg(1:97),umbcamp.WVfrac_percent_amp_sig,'k',0.25);
    hold off
    set(gca,'ydir','reverse'); plotaxis2;
    title('WVfrac\_percent amplification dWVfrac/dST'); ylim([100 1000]); xlim([-1 +1]*25); 
    hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);

  iFig = iFig + 1;
  figure(iFig); clf; 
  i300 = 1 : length(era5amp.T_amp);
  i300 = find(pavg >= 200); i300 = i300(1);
  i300 = find(pavg >= 200 & pavg <= 500); 
  if length(i300) == 1
    plot(era5amp.T_amp(i300),era5amp.WVfrac_percent_amp(i300),'bs',merra2amp.T_amp(i300),merra2amp.WVfrac_percent_amp(i300),'cs',...
         airsL3amp.T_amp(i300),airsL3amp.WVfrac_percent_amp(i300),'ro',climcapsL3amp.T_amp(i300),climcapsL3amp.WVfrac_percent_amp(i300),'mo',...
         umbcamp.T_amp(i300),umbcamp.WVfrac_percent_amp(i300),'k+','markersize',5,'linewidth',2)
    title('at 300 mb')
  else
    plot(nanmean(era5amp.T_amp(i300)),nanmean(era5amp.WVfrac_percent_amp(i300)),'b^',nanmean(merra2amp.T_amp(i300)),nanmean(merra2amp.WVfrac_percent_amp(i300)),'c^',...
         nanmean(airsL3amp.T_amp(i300)),nanmean(airsL3amp.WVfrac_percent_amp(i300)),'rx',nanmean(climcapsL3amp.T_amp(i300)),nanmean(climcapsL3amp.WVfrac_percent_amp(i300)),'mx',...
         nanmean(umbcamp.T_amp(i300)),nanmean(umbcamp.WVfrac_percent_amp(i300)),'kx','markersize',10,'linewidth',4)
    title('<200-500> mb')
  end
    xlabel('dT/dTs (K/K)'); ylabel('dq/dTs (%/K)'); plotaxis2;
    %ylim([0 20])
    %axis([0 2 0 20])
    hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10); 

  iFig = iFig + 1;
  figure(iFig); clf; 
    wonk = junkAIRSL3.ptemprate .* maskLF100;     miaow_T_AIRSL3      = 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkCLIMCAPSL3.ptemprate .* maskLF100; miaow_T_CLIMCAPSL3  = 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkERA5.ptemprate .* maskLF100;       miaow_T_ERA5        = 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkMERRA2.ptemprate .* maskLF100;     miaow_T_MERRA2      = 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkUMBC.ptemprate .* maskLF100;       miaow_T_UMBC        = 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkAIRSL3.waterrate .* maskLF100;     miaow_WV_AIRSL3     = 100 * 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkCLIMCAPSL3.waterrate .* maskLF100; miaow_WV_CLIMCAPSL3 = 100 * 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkERA5.waterrate .* maskLF100;       miaow_WV_ERA5       = 100 * 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkMERRA2.waterrate .* maskLF100;     miaow_WV_MERRA2     = 100 * 10 * nanmean(nanmean(wonk(i300,:),1));
    wonk = junkUMBC.waterrate .* maskLF100;       miaow_WV_UMBC       = 100 * 10 * nanmean(nanmean(wonk(i300,:),1));
    plot(miaow_T_ERA5,miaow_WV_ERA5,'b^',miaow_T_MERRA2,miaow_WV_MERRA2,'c^',...
         miaow_T_AIRSL3,miaow_WV_AIRSL3,'rx',miaow_T_CLIMCAPSL3,miaow_WV_CLIMCAPSL3,'mx',...
         miaow_T_UMBC,miaow_WV_UMBC,'kx','markersize',10,'linewidth',4)
    xlabel('dT/dt (K/decade)'); ylabel('dq/dt pc/decade)')
    hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS L3','UMBC','location','best','fontsize',10);


end
