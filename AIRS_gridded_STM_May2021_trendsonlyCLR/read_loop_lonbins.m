iDoAgain = +1;
iCountStart0 = 0;
iCount = 0;

while iDoAgain > 0
  for ii = 1 : 72
    iCountStart0 = iCount;  %% count at beginning of lonbin
    for jj = 1 : 64
      iibin = (jj-1)*72 + ii;
      if iOCBset == 0
        if dataset ~= 3
          if iNorD > 0
            fname = ['/asl/s1/sergio/Tiles4608/Output_WORKS_May18_2021_Great_AIRS_STM/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iibin) '.mat']; %% stored here after July 2021
            fname = ['Output/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iibin) '.mat']; %% stored here before before July 2021, and fornew test comparisons
          elseif iNorD < 0
            fname = ['Output_Day/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iibin) '.mat']; %% stored here before before July 2021
          end
        elseif dataset == 3
          if iNorD > 0
            fname = ['Output/Extreme/test' num2str(iibin) '.mat']; %% stored here before before July 2021, and fornew test comparisons
          elseif iNorD < 0
            fname = ['Output_Day/Extreme/test' num2str(iibin) '.mat']; %% stored here before before July 2021
          end
        elseif dataset == -3
          if iNorD > 0
            fname = ['Output/Quantile00/test' num2str(iibin) '.mat']; %% stored here before before July 2021, and fornew test comparisons
          elseif iNorD < 0
            fname = ['Output_Day/Quantile00/test' num2str(iibin) '.mat']; %% stored here before before July 2021
          end
        end
      elseif iOCBset == 1
        if iNorD > 0
          fname =     ['Output_CAL/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iibin) '.mat']; %% stored here before before July 2021, and fornew test comparisons
        elseif iNorD < 0
          fname = ['Output_Day_CAL/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iibin) '.mat']; %% stored here before before July 2021
        end
      end
      
      existfname(iibin) = exist(fname);
      i10sec = -1;
      if exist(fname)
        %% datenum : A serial date number represents the whole and fractional number of days from a fixed, preset date (January 0, 0000) in the proleptic ISO calendar.
        moo = dir(fname);
        rightnow = datenum(datetime('now'));
        if (rightnow-moo.datenum)*24*60*60 > 10
          i10sec = 1;
        end
      end
      
      %fprintf(1,'%5i    %s \n',existfname(iibin),fname)
      if exist(fname) > 0 & iaFound(iibin) == 0 & i10sec == 1
        fnamelastloaded = fname;
        iExist = +1;
        loader = ['load ' fname];
        eval(loader);
        iaFound(iibin) = +1;
        
        results(iibin,1:6)    = oem.finalrates(1:6);
        resultsunc(iibin,1:6) = oem.finalsigs(1:6);
        [mmn,nn] = size(oem.ak_water);
      
        if length(jacobian.wvjaclays_used) == iNumLay
          for iibini = 1 : length(jacobian.wvjaclays_used)
            %iavg = jacobian.wvjaclays_used{iNumLay}-6;
            iavg = jacobian.wvjaclays_used{iibini}-jacobian.wvjaclays_offset;
            pavg(iibini) = mean(plays(iavg));
          end
        end
      
        save_cov_set.cov_set(:,iibin)   = oem.cov_set;
        save_cov_set.fmat(:,iibin)      = sqrt(diag(oem.fmat));
        save_cov_set.reg_type        = oem.reg_type;
      
        %%% so six elements are set here : [1 2 3 4 5 6] typically from [CO2 N2O CH4 CFC11 CFC12 ST]
          save_cov_set.xb_trace(:,iibin)  = oem.xb(1:6);
          save_cov_set.xf_trace(:,iibin)  = oem.finalrates(1:6);
        %%% so two elements are set here : [1 2] typically from [1 length(jacobian.T/WV/O3_i)]
          save_cov_set.xb_traceI(:,iibin) = jacobian.scalar_i([1 length(jacobian.scalar_i)]); %% Co2/N2O/CH4/CFC11/CFC11/stemp etc
          save_cov_set.xb_wvzI(:,iibin)   = jacobian.water_i([1 length(jacobian.water_i)]);   %% top/bottom
          save_cov_set.xb_tzzI(:,iibin)   = jacobian.temp_i([1 length(jacobian.temp_i)]);     %% top/bottom
          save_cov_set.xb_ozzI(:,iibin)   = jacobian.ozone_i([1 length(jacobian.ozone_i)]);   %% top/bottom
          %
          save_cov_set.xb_wvz(:,iibin)    = oem.xb(jacobian.water_i([1 length(jacobian.water_i)]));  %% top/bottom
          save_cov_set.xb_tzz(:,iibin)    = oem.xb(jacobian.temp_i([1 length(jacobian.temp_i)]));    %% top/bottom
          save_cov_set.xb_ozz(:,iibin)    = oem.xb(jacobian.ozone_i([1 length(jacobian.ozone_i)]));  %% top/bottom
          %
          save_cov_set.xf_wvz(:,iibin)    = oem.finalrates(jacobian.water_i([1 length(jacobian.water_i)]));  %% top/bottom
          save_cov_set.xf_tzz(:,iibin)    = oem.finalrates(jacobian.temp_i([1 length(jacobian.temp_i)]));    %% top/bottom
          save_cov_set.xf_ozz(:,iibin)    = oem.finalrates(jacobian.ozone_i([1 length(jacobian.ozone_i)]));  %% top/bottom
      
        nlays_straight_from_results(iibin) = nn;
        nn0 = min(nn,iNumLay);
      
        xb(iibin,1:6)    = oem.xb(1:6);
        xbWV(iibin,1:nn) = oem.xb((1:nn)+6+nn*0);
        xbT(iibin,1:nn)  = oem.xb((1:nn)+6+nn*1);
        xbO3(iibin,1:nn) = oem.xb((1:nn)+6+nn*2);
      
        if nn0 == iNumLay
          thedofs(iibin) = oem.dofs;
          lencdofs(iibin) = length(oem.cdofs);
          cdofs(iibin,1:length(oem.cdofs)) = oem.cdofs;
          
          resultsWV(iibin,1:nn) = oem.finalrates((1:nn)+6+nn*0);
          resultsT(iibin,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
          resultsO3(iibin,1:nn) = oem.finalrates((1:nn)+6+nn*2);
          resultsWVunc(iibin,1:nn) = oem.finalsigs((1:nn)+6+nn*0);
          resultsTunc(iibin,1:nn)  = oem.finalsigs((1:nn)+6+nn*1);
          resultsO3unc(iibin,1:nn) = oem.finalsigs((1:nn)+6+nn*2);
      
        else
          iWarning = iWarning + 1;
          iaWarning(iWarning) = iibin;
      
          thedofs(iibin) = oem.dofs;
          lencdofs(iibin) = length(oem.cdofs);
          cdofs(iibin,1:length(oem.cdofs)) = oem.cdofs;
      
          resultsWV(iibin,:) = NaN;
          resultsT(iibin,:)  = NaN;
          resultsO3(iibin,:) = NaN;
          resultsWVunc(iibin,:) = NaN;
          resultsTunc(iibin,:)  = NaN;
          resultsO3unc(iibin,:) = NaN;
      
          wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(iibin,1:nn0) = wah(1:nn0);
          wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(iibin,1:nn0)  = wah(1:nn0);
          wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(iibin,1:nn0) = wah(1:nn0);
      
          wah = oem.finalsigs((1:nn)+6+nn*0);   resultsWVunc(iibin,1:nn0) = wah(1:nn0);
          wah = oem.finalsigs((1:nn)+6+nn*1);   resultsTunc(iibin,1:nn0)  = wah(1:nn0);
          wah = oem.finalsigs((1:nn)+6+nn*2);   resultsO3unc(iibin,1:nn0) = wah(1:nn0);
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
      
          %figure(47); plot(oem.ak_water',pjunk20,'c',max(oem.ak_water'),pjunk20,'rx-',mean(oem.ak_water'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
          %figure(48); plot(oem.ak_ozone',pjunk20,'c',max(oem.ak_ozone'),pjunk20,'rx-',mean(oem.ak_ozone'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
          %figure(49); plot(oem.ak_temp',pjunk20,'c',max(oem.ak_temp'),pjunk20,'rx-',mean(oem.ak_temp'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
      
          clear waterrate_ak0 waterrate_ak1 o3rate_ak0 o3rate_ak1 temprate_ak0 temprate_ak1
          ix = iibin;
          waterrate_ak0 = ones(ix,1)*era5.trend_gas_1(1:100,ix)';
            for iibini = 1 : length(jacobian.wvjaclays_used)
              junk = jacobian.wvjaclays_used{iibini}-6;
              waterrate_ak1(:,iibini) = mean(waterrate_ak0(:,junk)');
            end
            ak = oem.ak_water;
            mean_ak_wv(ix,1:length(ak)) = max(ak');
            waterrate_ak0_era5(ix,1:length(ak)) = (waterrate_ak1(ix,:)')';
            waterrate_akF_era5(ix,1:length(ak)) = (ak * waterrate_ak1(ix,:)')';
          o3rate_ak0 = ones(ix,1)*era5.trend_gas_3(1:100,ix)';
            for iibini = 1 : length(jacobian.wvjaclays_used)
              junk = jacobian.wvjaclays_used{iibini}-6;
              o3rate_ak1(:,iibini) = mean(o3rate_ak0(:,junk)');
            end
            ak = oem.ak_ozone;
            mean_ak_o3(ix,1:length(ak)) = max(ak');
            o3rate_ak0_era5(ix,1:length(ak)) = (o3rate_ak1(ix,:)')';
            o3rate_akF_era5(ix,1:length(ak)) = (ak * o3rate_ak1(ix,:)')';
          temprate_ak0 = ones(ix,1)*era5.trend_ptemp(1:100,ix)';
            for iibini = 1 : length(jacobian.wvjaclays_used)
              junk = jacobian.wvjaclays_used{iibini}-6;
              temprate_ak1(:,iibini) = mean(temprate_ak0(:,junk)');
            end
            ak = oem.ak_temp;
            mean_ak_T(ix,1:length(ak)) = max(ak');
            temprate_ak0_era5(ix,1:length(ak)) = (temprate_ak1(ix,:)')';
            temprate_akF_era5(ix,1:length(ak)) = (ak * temprate_ak1(ix,:)')';
        end   %% if iAK > 0
      
        %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%%
      
        junknoise  = nan(iNumChan,1);
        junknoise2 = nan(iNumChan,1);
        junknoise(jacobian.chanset)  = sqrt(diag(oem.se));
        junknoise2(jacobian.chanset) = rateset.unc_rates(jacobian.chanset);
        junknoise2                   = rateset.unc_rates;
      
        %% depending on the lag1 correction, maybe junknoise and junknoise2 are different
        if sum([junknoise(jacobian.chanset)-junknoise2(jacobian.chanset)]) > eps
          fprintf(1,'WARNING fov %4i we have junknoise and junknoise2 differing \n',iibin)
        end
      
        if isfield(oem,'spectral_deltan00')
          spectral_deltan00(:,iibin) = oem.spectral_deltan00;
        end
      
        rates(:,iibin)           = rateset.rates;
        fits(:,iibin)            = oem.fit';
        if isfield(oem,'fitXcomponents')
          componentfits(1,:,iibin) = oem.fitXcomponents(1,:);   %% trace gases CO2/N2O/CH4
          componentfits(2,:,iibin) = oem.fitXcomponents(2,:);   %% ST
          componentfits(3,:,iibin) = oem.fitXcomponents(3,:);   %% WV(z)
          componentfits(4,:,iibin) = oem.fitXcomponents(4,:);   %% T(z)
          componentfits(5,:,iibin) = oem.fitXcomponents(5,:);   %% O3(z)
        end
      
        nedt(:,iibin)  = junknoise;  %% will have lots of NaNs since based on oem.se which only used selected channels
        nedt(:,iibin)  = junknoise2; %% should be nicely filled in
      
      elseif exist(fname) > 0 & iaFound(iibin) == 1
        %% do nothing, things are cool
      else
        iExist = -1;
        iaFound(iibin) = 0;
        results(iibin,1:6) = NaN;
        resultsWV(iibin,:) = NaN;
        resultsT(iibin,:)  = NaN;
        resultsO3(iibin,:) = NaN;
      
        resultsunc(iibin,1:6) = NaN;
        resultsWVunc(iibin,:) = NaN;
        resultsTunc(iibin,:)  = NaN;
        resultsO3unc(iibin,:) = NaN;
      
        rates(:,iibin) = NaN;
        fits(:,iibin) = NaN;
        nedt(:,iibin) = NaN;
      end
    end    %% for jj loop over lats 1--64

    if mod(ii,10) == 0
      fprintf(1,'+')
    else
      fprintf(1,'.');
    end
    fprintf(1,'lonbin = %2i ... looped over all 64 latbins \n',ii)
    
    figure(6)
    pcolor(reshape(results(:,6),72,64)); colorbar; caxis([-1 +1]*0.15); colormap(llsmap5);

    iCount = sum(iaFound);
    if iCount > iCountStart0
      %% do the plots
      load latB64.mat
      rlat65 = latB2; rlon73 = -180 : 5 : +180;
      rlon = -180 : 5 : +180;  rlat = latB2; 
      rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
      rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
      aslmap(6,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt so far');     caxis([-1 +1]*0.15); colormap(llsmap5)
      
      jett = jet(64); jett(1,:) = 1;
      %figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(rlat,1:iNumLay,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dT/dt');      colormap(llsmap5); caxis([-1 +1]*0.15)
      %figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(rlat,1:iNumLay,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.015)
      figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(rlat,pavg,waha);  shading interp; colorbar('horizontal'); set(gca,'ydir','reverse'); title('UMBC dT/dt');      colormap(llsmap5); caxis([-1 +1]*0.15)
      figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(rlat,pavg,waha);  shading interp; colorbar('horizontal'); set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.01)
        figure(29); set(gca,'yscale','log'); ylim([10 1000]);       figure(30); set(gca,'yscale','linear'); ylim([100 1000]);
      figure(31); clf; waha = reshape(iaFound,72,64);                                                   pcolor(rlon,rlat,waha'); shading flat;   colorbar; set(gca,'ydir','normal');  
        title([num2str(sum(iaFound(:))) ' / 4608 = ' num2str(100*sum(iaFound(:))/4608) ' % made so far']);  
        xlabel('Longitude'); ylabel('Latitude'); colormap(jett); 
      pause(0.1)
    end
    
  end      %% for ii loop over lons 1-72

  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  fprintf(1,'\n');
  resultsunc = real(resultsunc);
  resultsWVunc = real(resultsWVunc);
  resultsTunc  = real(resultsTunc);
  resultsO3unc = real(resultsO3unc);
  
  display('Trying to look at atributes of the last file that should have been read in .... ')
  lser = ['!ls -lt ' fnamelastloaded]; eval(lser)

  iCount = sum(iaFound);
  fprintf(1,'found %4i of %4i \n',sum(iaFound),64*72)
  iDoAgain = -1;
  if sum(iaFound) < 64*72
    figure(6)
    pcolor(reshape(results(:,6),72,64)); colorbar; caxis([-1 +1]*0.15); colormap(llsmap5);

    load latB64.mat
    rlat65 = latB2; rlon73 = -180 : 5 : +180;
    rlon = -180 : 5 : +180;  rlat = latB2; 
    rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
    rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
    aslmap(6,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt so far');     caxis([-1 +1]*0.15); colormap(llsmap5)
    
    jett = jet(64); jett(1,:) = 1;
    %figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(rlat,1:iNumLay,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dT/dt');      colormap(llsmap5); caxis([-1 +1]*0.15)
    %figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(rlat,1:iNumLay,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.015)
    figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(rlat,pavg,waha);  shading interp; colorbar('horizontal'); set(gca,'ydir','reverse'); title('UMBC dT/dt');      colormap(llsmap5); caxis([-1 +1]*0.15)
    figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(rlat,pavg,waha);  shading interp; colorbar('horizontal'); set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.01)
      figure(29); set(gca,'yscale','log'); ylim([10 1000]);       figure(30); set(gca,'yscale','linear'); ylim([100 1000]);
    figure(31); clf; waha = reshape(iaFound,72,64);                                                   pcolor(rlon,rlat,waha'); shading flat;   colorbar; set(gca,'ydir','normal');  
      title([num2str(sum(iaFound(:))) ' / 4608 = ' num2str(100*sum(iaFound(:))/4608) ' % made so far']);  
      xlabel('Longitude'); ylabel('Latitude'); colormap(jett); 
    iDoAgain = input('read in remaining files (-1/+1 Default) : '); 
    if length(iDoAgain) == 0
      iDoAgain = +1;
    end
  end
end
