fprintf(1,'will loop through %3i timeteps for %2i latbins .. the timesteps will mark off as "o" for 100 and "." for 10 \n',iNumAnomTimeSteps,iNumAnomTiles)
iaaFound = nan(iNumAnomTimeSteps,iNumAnomTiles);
while iDoAgain > 0
  for iix = 1 : iNumAnomTimeSteps    
    for iiy = 1 : iNumAnomTiles
      ii = (iiy-1)*iNumAnomTimeSteps + iix;
      iInd = ii;
      set_anom_outfilename        
      
      fname = [zanom_outdir '/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    
      existfname(ii) = exist(fname);
      i10sec = -1;
      if exist(fname)
        %% datenum : A serial date number represents the whole and fractional number of days from a fixed, preset date (January 0, 0000) in the proleptic ISO calendar.
        moo = dir(fname);
        rightnow = datenum(datetime('now'));
        if (rightnow-moo.datenum)*24*60*60 > 10
          i10sec = 1;
        end
      end
      %fprintf(1,'%5i    %s \n',existfname(ii),fname)
      if exist(fname) > 0 & iaFound(ii) == 0 & i10sec == 1
        fnamelastloaded = fname;
        iExist = +1;
        loader = ['load ' fname];
        eval(loader);
	iaaFound(iix,iiy) = 1;
        iaFound(ii) = +1;
        
        results(ii,1:6)    = oem.finalrates(1:6);
        resultsunc(ii,1:6) = oem.finalsigs(1:6);
        [mmn,nn] = size(oem.ak_water);
  
        if length(jacobian.wvjaclays_used) == iNumLay
          for iii = 1 : length(jacobian.wvjaclays_used)
            %iavg = jacobian.wvjaclays_used{iNumLay}-6;
            iavg = jacobian.wvjaclays_used{iii}-jacobian.wvjaclays_offset;
            pavg(iii) = mean(plays(iavg));
          end
        end
  
        save_cov_set.cov_set(:,ii)   = oem.cov_set;
        save_cov_set.fmat(:,ii)      = sqrt(diag(oem.fmat));
        save_cov_set.reg_type        = oem.reg_type;
        save_cov_set.xb_traceI(:,ii) = jacobian.scalar_i([1 length(jacobian.scalar_i)]); %% Co2/N2O/CH4/CFC11/CFC11/stemp etc  indices (1,6)
        save_cov_set.xb_wvzI(:,ii)   = jacobian.water_i([1 length(jacobian.water_i)]);   %% top/bottom indices
        save_cov_set.xb_tzzI(:,ii)   = jacobian.temp_i([1 length(jacobian.temp_i)]);     %% top/bottom indices
        save_cov_set.xb_ozzI(:,ii)   = jacobian.ozone_i([1 length(jacobian.ozone_i)]);   %% top/bottom indices
        save_cov_set.xb_trace(:,ii)  = oem.xb(1:6);                                             %% Co2/N2O/CH4/CFC11/CFC11/stemp etc  init values (1,6)
        save_cov_set.xb_wvz(:,ii)    = oem.xb(jacobian.water_i([1 length(jacobian.water_i)]));  %% top/bottom init values
        save_cov_set.xb_tzz(:,ii)    = oem.xb(jacobian.temp_i([1 length(jacobian.temp_i)]));    %% top/bottom init values
        save_cov_set.xb_ozz(:,ii)    = oem.xb(jacobian.ozone_i([1 length(jacobian.ozone_i)]));  %% top/bottom init values
        save_cov_set.xf_trace(:,ii)  = oem.finalrates(1:6);                                             %% Co2/N2O/CH4/CFC11/CFC11/stemp etc  final values (1,6)
        save_cov_set.xf_wvz(:,ii)    = oem.finalrates(jacobian.water_i([1 length(jacobian.water_i)]));  %% top/bottom final values
        save_cov_set.xf_tzz(:,ii)    = oem.finalrates(jacobian.temp_i([1 length(jacobian.temp_i)]));    %% top/bottom final values
        save_cov_set.xf_ozz(:,ii)    = oem.finalrates(jacobian.ozone_i([1 length(jacobian.ozone_i)]));  %% top/bottom final values
  
        nlays_straight_from_results(ii) = nn;
        nn0 = min(nn,iNumLay);
  
        rtime(ii)     = anomalyinfo.rtime;
        xb(ii,1:6)    = oem.xb(1:6);
        xbWV(ii,1:nn) = oem.xb((1:nn)+6+nn*0);
        xbT(ii,1:nn)  = oem.xb((1:nn)+6+nn*1);
        xbO3(ii,1:nn) = oem.xb((1:nn)+6+nn*2);
  
        if nn0 == iNumLay
          thedofs(ii) = oem.dofs;
          lencdofs(ii) = length(oem.cdofs);
          cdofs(ii,1:length(oem.cdofs)) = oem.cdofs;
          
          resultsWV(ii,1:nn) = oem.finalrates((1:nn)+6+nn*0);
          resultsT(ii,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
          resultsO3(ii,1:nn) = oem.finalrates((1:nn)+6+nn*2);
          resultsWVunc(ii,1:nn) = oem.finalsigs((1:nn)+6+nn*0);
          resultsTunc(ii,1:nn)  = oem.finalsigs((1:nn)+6+nn*1);
          resultsO3unc(ii,1:nn) = oem.finalsigs((1:nn)+6+nn*2);
  
        else
          iWarning = iWarning + 1;
          iaWarning(iWarning) = ii;
    
          thedofs(ii) = oem.dofs;
          lencdofs(ii) = length(oem.cdofs);
          cdofs(ii,1:length(oem.cdofs)) = oem.cdofs;
    
          resultsWV(ii,:) = NaN;
          resultsT(ii,:)  = NaN;
          resultsO3(ii,:) = NaN;
          resultsWVunc(ii,:) = NaN;
          resultsTunc(ii,:)  = NaN;
          resultsO3unc(ii,:) = NaN;
    
          wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(ii,1:nn0) = wah(1:nn0);
          wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(ii,1:nn0)  = wah(1:nn0);
          wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(ii,1:nn0) = wah(1:nn0);
    
          wah = oem.finalsigs((1:nn)+6+nn*0);   resultsWVunc(ii,1:nn0) = wah(1:nn0);
          wah = oem.finalsigs((1:nn)+6+nn*1);   resultsTunc(ii,1:nn0)  = wah(1:nn0);
          wah = oem.finalsigs((1:nn)+6+nn*2);   resultsO3unc(ii,1:nn0) = wah(1:nn0);
        end
  
        %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%%
        if iAK > 0
          if ~exist('waterrate_akF_era5')
             waterrate_akF_era5 = nan(size(resultsWV));
             o3rate_akF_era5 = nan(size(resultsWV));
             temprrate_akF_era5 = nan(size(resultsWV));
  
             waterrate_ak0_era5 = nan(size(resultsWV));
             o3rate_ak0_era5 = nan(size(resultsWV));
             temprrate_ak0_era5 = nan(size(resultsWV));
  
             mean_ak_wv = nan(size(resultsWV));;
             mean_ak_T  = nan(size(resultsWV));;
             mean_ak_o3 = nan(size(resultsWV));;
          end
    
          %figure(2); plot(oem.ak_water',pjunk20,'c',max(oem.ak_water'),pjunk20,'rx-',mean(oem.ak_water'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
          %figure(3); plot(oem.ak_ozone',pjunk20,'c',max(oem.ak_ozone'),pjunk20,'rx-',mean(oem.ak_ozone'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
          %figure(4); plot(oem.ak_temp',pjunk20,'c',max(oem.ak_temp'),pjunk20,'rx-',mean(oem.ak_temp'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
  
          clear waterrate_ak0 waterrate_ak1 o3rate_ak0 o3rate_ak1 temprate_ak0 temprate_ak1
          ix = ii;
          waterrate_ak0 = ones(ix,1)*era5.trend_gas_1(1:100,ix)';
            for iii = 1 : length(jacobian.wvjaclays_used)
              junk = jacobian.wvjaclays_used{iii}-6;
              waterrate_ak1(:,iii) = mean(waterrate_ak0(:,junk)');
            end
            ak = oem.ak_water;
            mean_ak_wv(ix,1:length(ak)) = max(ak');
            waterrate_ak0_era5(ix,1:length(ak)) = (waterrate_ak1(ix,:)')';
            waterrate_akF_era5(ix,1:length(ak)) = (ak * waterrate_ak1(ix,:)')';
          o3rate_ak0 = ones(ix,1)*era5.trend_gas_3(1:100,ix)';
            for iii = 1 : length(jacobian.wvjaclays_used)
              junk = jacobian.wvjaclays_used{iii}-6;
              o3rate_ak1(:,iii) = mean(o3rate_ak0(:,junk)');
            end
            ak = oem.ak_ozone;
            mean_ak_o3(ix,1:length(ak)) = max(ak');
            o3rate_ak0_era5(ix,1:length(ak)) = (o3rate_ak1(ix,:)')';
            o3rate_akF_era5(ix,1:length(ak)) = (ak * o3rate_ak1(ix,:)')';
          temprate_ak0 = ones(ix,1)*era5.trend_ptemp(1:100,ix)';
            for iii = 1 : length(jacobian.wvjaclays_used)
              junk = jacobian.wvjaclays_used{iii}-6;
              temprate_ak1(:,iii) = mean(temprate_ak0(:,junk)');
            end
            ak = oem.ak_temp;
            mean_ak_T(ix,1:length(ak)) = max(ak');
            temprate_ak0_era5(ix,1:length(ak)) = (temprate_ak1(ix,:)')';
            temprate_akF_era5(ix,1:length(ak)) = (ak * temprate_ak1(ix,:)')';
        end   %% if iAK > 0
        %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%%
  
        junknoise  = nan(2645,1);
        junknoise2 = nan(2645,1);
        junknoise(jacobian.chanset)  = sqrt(diag(oem.se));
        junknoise2(jacobian.chanset) = rateset.unc_rates(jacobian.chanset);
        junknoise2                   = rateset.unc_rates;
  
        %% depending on the lag1 correction, maybe junknoise and junknoise2 are different
        if sum([junknoise(jacobian.chanset)-junknoise2(jacobian.chanset)]) > eps
          fprintf(1,'WARNING fov %4i we have junknoise and junknoise2 differing \n',ii)
        end
    
        if isfield(oem,'spectral_deltan00')
          spectral_deltan00(:,ii) = oem.spectral_deltan00;
        end
  
        rates(:,ii)           = rateset.rates;
        fits(:,ii)            = oem.fit';
        if isfield(oem,'fitXcomponents')
          componentfits(1,:,ii) = oem.fitXcomponents(1,:);   %% trace gases CO2/N2O/CH4
          componentfits(2,:,ii) = oem.fitXcomponents(2,:);   %% ST
          componentfits(3,:,ii) = oem.fitXcomponents(3,:);   %% WV(z)
          componentfits(4,:,ii) = oem.fitXcomponents(4,:);   %% T(z)
          componentfits(5,:,ii) = oem.fitXcomponents(5,:);   %% O3(z)
        end
  
        nedt(:,ii)  = junknoise;  %% will have lots of NaNs since based on oem.se which only used selected channels
        nedt(:,ii)  = junknoise2; %% should be nicely filled in
  
      elseif exist(fname) > 0 & iaFound(ii) == 1
        %% do nothing, things are cool
      else
        iExist = -1;
        iaFound(ii) = 0;
        results(ii,1:6) = NaN;
        resultsWV(ii,:) = NaN;
        resultsT(ii,:)  = NaN;
        resultsO3(ii,:) = NaN;
    
        resultsunc(ii,1:6) = NaN;
        resultsWVunc(ii,:) = NaN;
        resultsTunc(ii,:)  = NaN;
        resultsO3unc(ii,:) = NaN;
    
        rates(:,ii) = NaN;
        fits(:,ii) = NaN;
        nedt(:,ii) = NaN;
      end

    end      %% loop over latbins
    figure(1);
      pcolor(iaaFound'); shading interp; colorbar; colormap jet;
      title(['iaaFound = ' num2str(sum(iaFound)) '/' num2str(iNumAnomTimeSteps*iNumAnomTiles) ' = ' num2str(100*sum(iaFound)/(iNumAnomTimeSteps*iNumAnomTiles)) ' %']);
    figure(2); pcolor(reshape(results(:,6),iNumAnomTimeSteps,iNumAnomTiles)'); title('SKT anomaly'); shading interp; colorbar; colormap(usa2); caxis([-1 +1]/250);
      title('SKT anomaly')
    pause(0.1)
    
  end    %% loop over timesteps
  
  fprintf(1,'\n');
  resultsunc = real(resultsunc);
  resultsWVunc = real(resultsWVunc);
  resultsTunc  = real(resultsTunc);
  resultsO3unc = real(resultsO3unc);
  
  display('Trying to look at atributes of the last file that should have been read in .... ')
  lser = ['!ls -lt ' fnamelastloaded]; eval(lser)
  
  fprintf(1,'found %4i of %4i \n',sum(iaFound),iNumAnomData)
  iDoAgain = -1;

  if sum(iaFound) < iNumAnomData
    if strfind(anomalydatafile,'_tile')
      plot_anomalies_1
    else
      simple = +1;
      plot_anomalies_All
    end

    junk_globalavg_rawdata = data_anom.btavgAnomFinal(:,1:iNumAnomTimeSteps);
    junk_globalavg_spectral_deltan00 = spectral_deltan00(:,1:iNumAnomTimeSteps);
    figure(3); clf; 
    plot(yymm,smooth(junk_globalavg_spectral_deltan00(i0723,:),23),'r.-',...
        yymm,smooth(junk_globalavg_rawdata(i0723,:),23),'r--',...
       'linewidth',2); plotaxis2;
    title('Fitted BT0723 data in thick \newline raw data in dashes, CO2 only'); pause(0.1)

    figure(4); clf
    hold off; wah = save_cov_set.xf_trace(1:6,1:iNumAnomTimeSteps);
      plot(yymm,wah(1,:),'b',yymm,wah(2,:),'m',yymm,wah(3,:),'g',yymm,wah(4,:),'y',yymm,wah(5,:),'k',yymm,wah(6,:),'r','linewidth',2); hold on; 
    hold on; wah = save_cov_set.xb_trace(1:6,1:iNumAnomTimeSteps);
      plot(yymm,wah(1,:),'b--',yymm,wah(2,:),'m--',yymm,wah(3,:),'g--',yymm,wah(4,:),'y--',yymm,wah(5,:),'k--',yymm,wah(6,:),'r--','linewidth',2); hold on; 
      hl = legend('CO2','N2O','CH4','CFC11','CFC12','ST','location','best');
      set(gca,'fontsize',10); title('xb(trace gas + ST)');
    hold off

    iDoAgain = input('read in remaining files (-1/+1 Default) : '); 
    if length(iDoAgain) == 0
      iDoAgain = +1;
    end
  end
end
