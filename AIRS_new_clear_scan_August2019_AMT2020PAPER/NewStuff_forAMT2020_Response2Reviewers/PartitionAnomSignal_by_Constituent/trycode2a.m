%% based on driver_plot_T_WV_anomaly_retrievals.m

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP

latbins = equal_area_spherical_bands(21);
latbinsx = 0.5*(latbins(1:end-1)+latbins(2:end));
tropics = find(abs(latbinsx) <= 30);
tropics = 11 : 30;
extralats_tropics = find(abs(latbinsx) <= 50);   %% this is 6 - 37

iOBSorCAL = 0;
iAvgNumDays = 0;
dataset = 1; %% AIRS

iLatbin = 20;
if iOBSorCAL == 0 & dataset == 1
  fin = ['../../ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 1
  fin = ['../../ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
elseif iOBSorCAL == 0 & dataset == 2
  fin = ['../../IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 2
  fin = ['../../IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
end

load(fin)
save_days = avg_doy_since2002;
save_rtime = avg16_rtime;
save_dat_1231 = zeros(length(save_days),40);
save_days_map = zeros(length(save_days),40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file1 = 'NewStuff_forAMT2020_Response2Reviewers/PartitionAnomSignal_by_Constituent/complete_jac_anomaly_0dayavg_results';
get_20lays = load(file1);
disp('loading up the tropical data bin by bin, timestep by timestep')
if ~exist('jacmeanTrop','var') & ~exist('ind_components','var')
  if ~exist('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/only_tropics_jacs_ind_components.mat')
    disp('doing trops only')
    iiStart = 1; iiEnd = 365;
    iiStart = 200; iiEnd = 200;
    iiStart = 300; iiEnd = 300;
    iiStart = 200; iiEnd = 300;
    iiStart = 001; iiEnd = 365;
    for ii = iiStart : iiEnd
      fprintf(1,'ii = %3i \n',ii)
      fname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' num2str(ii,'%03d') '_M_TS_jac_all_5_97_97_97_2645.mat'];
      load(fname)
      jacjunk = M_TS_jac_all(tropics,:,:);
      tracegas_jacjunk = jacjunk(:,:,1:6);
      for kk = 1 : 6
        tracegas_jacjunk2(:,:,kk) = tracegas_jacjunk(:,:,kk)/qrenorm(kk);
      end
      wv_jacjunk = jacjunk(:,:,(1:97)+0*97+6)/qrenorm(020);
      tz_jacjunk = jacjunk(:,:,(1:97)+1*97+6)/qrenorm(120);
      o3_jacjunk = jacjunk(:,:,(1:97)+2*97+6)/qrenorm(220);

      wv_jacjunk2 = zeros(length(tropics),2645,20);
      tz_jacjunk2 = zeros(length(tropics),2645,20);
      o3_jacjunk2 = zeros(length(tropics),2645,20);
      for kk = 1 : length(tropics)
        junk = get_20lays.a.jacobian.wvjaclays_used{kk}-6;
        wv_jacjunk2(:,:,kk) = sum(wv_jacjunk(:,:,junk),3);
        tz_jacjunk2(:,:,kk) = sum(tz_jacjunk(:,:,junk),3);
        o3_jacjunk2(:,:,kk) = sum(o3_jacjunk(:,:,junk),3);
      end

      for kk = 1 : length(tropics)
        xresults = load(['../../SAVE_LW_noCFC11_noN2O/OutputAnomaly_OBS/' num2str(tropics(kk),'%02d') '/anomtest_timestep' num2str(ii) '.mat']);
        ind_componentsTGST(kk,ii-iiStart+1,:,1) = xresults.oem.finalrates(1) * squeeze(tracegas_jacjunk2(kk,:,1));
        ind_componentsTGST(kk,ii-iiStart+1,:,2) = xresults.oem.finalrates(2) * squeeze(tracegas_jacjunk2(kk,:,2));
        ind_componentsTGST(kk,ii-iiStart+1,:,3) = xresults.oem.finalrates(3) * squeeze(tracegas_jacjunk2(kk,:,3));
        ind_componentsTGST(kk,ii-iiStart+1,:,4) = xresults.oem.finalrates(4) * squeeze(tracegas_jacjunk2(kk,:,4));
        ind_componentsTGST(kk,ii-iiStart+1,:,5) = xresults.oem.finalrates(5) * squeeze(tracegas_jacjunk2(kk,:,5));
        ind_componentsTGST(kk,ii-iiStart+1,:,6) = xresults.oem.finalrates(6) * squeeze(tracegas_jacjunk2(kk,:,6));

        junk = xresults.oem.finalrates(6+(1:20)+0*20);
        junk = (ones(2645,1) * junk') .* squeeze(wv_jacjunk2(kk,:,:));
        ind_componentsWV(kk,ii-iiStart+1,:) = sum(junk,2);

        junk = xresults.oem.finalrates(6+(1:20)+1*20);
        junk = (ones(2645,1) * junk') .* squeeze(tz_jacjunk2(kk,:,:));
        ind_componentsTZ(kk,ii-iiStart+1,:) = sum(junk,2);

        junk = xresults.oem.finalrates(6+(1:20)+2*20);
        junk = (ones(2645,1) * junk') .* squeeze(o3_jacjunk2(kk,:,:));
        ind_componentsO3(kk,ii-iiStart+1,:) = sum(junk,2);
        
        ind_raaObs(kk,ii-iiStart+1,:) = xresults.rateset.rates;
        ind_raaCal(kk,ii-iiStart+1,:) = xresults.oem.fit;

        ind_raaReconstruct(kk,ii-iiStart+1,:) = sum(ind_componentsTGST(kk,ii-iiStart+1,:,:),4) + ind_componentsWV(kk,ii-iiStart+1,:) + ...
                                                 ind_componentsTZ(kk,ii-iiStart+1,:) + ind_componentsO3(kk,ii-iiStart+1,:);
         plot(1:2645,squeeze(ind_raaObs(kk,ii-iiStart+1,:)),'r',1:2645,squeeze(ind_raaCal(kk,ii-iiStart+1,:)),'bo-',...
              1:2645,squeeze(ind_raaReconstruct(kk,ii-iiStart+1,:)),'k')
      end
    end
    save -v7.3 /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/only_tropics_ind_componentsA.mat f qrenorm ind_*
  else
    load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/only_tropics_ind_componentsA.mat
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raaObs = squeeze(nanmean(ind_raaObs,1));
raaCal = squeeze(nanmean(ind_raaCal,1));
raaRec = squeeze(nanmean(ind_raaReconstruct,1));

wv = squeeze(nanmean(ind_componentsWV,1));
tz = squeeze(nanmean(ind_componentsTZ,1));
o3 = squeeze(nanmean(ind_componentsO3,1));
tgst = squeeze(nanmean(ind_componentsTGST,1));

if (iiEnd > iiStart)
  raaObs = nanmean(raaObs,1);
  raaCal = nanmean(raaCal,1);
  raaRec = nanmean(raaRec,1);

  wv = nanmean(wv,1);
  tz = nanmean(tz,1);
  o3 = nanmean(o3,1);
  tgst = squeeze(nanmean(tgst,1));

end
figure(1)

plot(f,raaCal - raaRec,'g'); %% remember above I read in kCARTA analytic CO2 jacs, but we really did fiite diff CO2 jacs
plot(f,raaObs,'k','linewidth',2); hold; plot(f,raaCal,'b.-',f,raaRec,'c'); hold off

plot(f,raaObs,'k','linewidth',2); hold; plot(f,tgst(:,[1 2 3 6]),f,wv,f,tz,f,o3); hold off
hl = legend('signal','co2','n2o','ch4','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raaObs = squeeze(nanmean(ind_raaObs,1));
raaCal = squeeze(nanstd(ind_raaCal,1));
raaRec = squeeze(nanstd(ind_raaReconstruct,1));

wv = squeeze(nanmean(ind_componentsWV,1));
tz = squeeze(nanmean(ind_componentsTZ,1));
o3 = squeeze(nanmean(ind_componentsO3,1));
tgst = squeeze(nanmean(ind_componentsTGST,1));

if (iiEnd > iiStart)
  raaObs = nanstd(raaObs,1);
  raaCal = nanstd(raaCal,1);
  raaRec = nanstd(raaRec,1);

  wv = nanstd(wv,1);
  tz = nanstd(tz,1);
  o3 = nanstd(o3,1);
  tgst = squeeze(nanstd(tgst,1));

end

figure(2)
plot(f,raaCal - raaRec,'g'); %% remember above I read in kCARTA analytic CO2 jacs, but we really did fiite diff CO2 jacs
plot(f,raaObs,'k','linewidth',2); hold; plot(f,raaCal,'b.-',f,raaRec,'c'); hold off

plot(f,raaObs,'k','linewidth',2); hold; plot(f,tgst(:,[1 2 3 6]),f,wv,f,tz,f,o3); hold off
hl = legend('signal','co2','n2o','ch4','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
xlim([640 1640])
