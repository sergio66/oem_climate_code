settings.iaSequential = -1;                  %% default is to retrieve everything at one go, 
                                             %% or can do sequential ala Chris Barnet eg do [150 60 100] to do T WV O3, after "fixing" CO2 and CH4 and N2O

settings.iSergioCO2 = -1;                    %% assume ESRL CO2/CH4 rates, DEFAULT   (+1) to fit for CO2/CH4

settings.set_era5_cmip6_airsL3         = -1; %% default, no a priori, else set to 3,5,6 for AIRS L3/ERA5/CMIP6
settings.set_era5_cmip6_airsL3_WV_T_O3 = -1; %% is using ERA5 or AIRS L3 or MERRA to set rates, you can choose to see -1:WV/T/ST/O3 or +1/+2/+3/+4/+5/+40  for WV/T+ST/O3/T/ST/LowerT only

settings.resetnorm2one = -1; %%% default, keep my scaling, set to +1 if you want to reset all to 1.00000000

settings.iNoiseType = -1;   % use squeeze(b_err_desc(ix,iy,:)) since these are TRENDS
                            % +1 uses the ideas of AIRS stability paper = NeDT/sqrt(N) but those are aomalies

settings.dataset = 1;       % see driver_put_together_QuantileChoose_trends.m
                            % (+1) AIRS 17 year quantile dataset, Strow  March 2021 2002/09-2019/08
                            % (-1) AIRS 17 year quantile dataset, Sergio Aug   2021 2002/09-2019/08
                            % (+2) AIRS 19 year quantile dataset, Sergio Aug   2021 2002/09-2021/08
                            % (+3) AIRS 19 year extreme  dataset, Sergio Aug   2021 2002/09-2021/08
                            % (-3) AIRS 19 year mean     dataset, Sergio Aug   2021 2002/09-2021/08
                            
settings.co2lays = 1;       % assume column jac for CO2, or 3 lays (gnd-500,500-trop,trop-TOA

settings.ocb_set = 0;       % 0 = obs, 1 = cal, -1 = bias
settings.numchan = 2645;    % L1b = 2378, L1c = 2645
settings.chan_LW_SW = 0;    % 0 is 640 to 1640, 1 is 700 to 1640, -1 is 640 to 2740

settings.set_tracegas = -1;            %% do we leave apriori as 0 or set CO2/N2o/CH4/CFC to be 2.2, 4.5, -1             
                                       %%   (-1 = no, +1 = yes)
                                       %%   note : if anomaly, adjust a priori for CO2,N2O,CH4,CFC
settings.offsetrates  = -1;            %% do we add a constant offset to the spectral rates (-1 = no, +1 = yes)
                                       %% if driver.rateset.ocb_set  == 'obs';
settings.addco2jacs   = -1;            %% do we add co2 jacs to the spectral rates (-1 = no, +1 = yes)
                                       %% if driver.rateset.ocb_set  == 'cal';
settings.obs_corr_matrix = -1;         %% just use nc_error (-1) or try to be fancy and use full cov matrix (+1)

settings.tie_sst_lowestlayer = +1;     %% tie together SST with lowest T(z)
settings.invtype         = 1;          %% pinv, see /home/sergio/MATLABCODE/oem_pkg/rodgers.m
settings.iNlays_retrieve = 97;         %% do all 97 layers
settings.descORasc = +1;               %% descending default
settings.iXJac = 0;                    %% const geo jacs, replace as needed CO2/CH4/N20;  
                                       %% +2 uses kCARTA varying geo/trace, +1 uses SARTA varying geo/trace jacs, 0 = constant kcarta jacs
settings.iDoStrowFiniteJac = -1;       %% -1 : do not change the time varying anomaly jacs                                done for all anomaly timesteps
                                       %% +1 stick to Sergio tracegas jacs = BT(1.001 X(t,latbin)) - BT(1.00 X(t,latbin)) interp in time
                                       %% +2 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      interp in time .. 
                                       %% +3 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      done for all anomaly timesteps DEFAULT
                                       %% +4 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      done for all anomaly timesteps with Age of Air for CO2
settings.iChSet = 1;                   %% +1 default, old chans (about 500)
                                       %% +2, new chans (about 400) with CFC11,CFC12      and weak WV, bad chans gone
                                       %% +3, new chans (about 400) w/o  CFC11 with CFC12 and weak WV, bad chans gone

settings.iFixTz_NoFit = -1;            %% -1 : do not fix Tz to ERA anomaly T(z,time) values, then fit for Tz
                                       %% +1 : do     fix Tz to ERA anomaly T(z,time) values, then keep Tz fixed (ie do not fit)
settings.iFixWV_NoFit = -1;            %% -1 : do not fix WV to ERA anomaly T(z,time) values, then fit for WV
                                       %% +1 : do     fix WV to ERA anomaly T(z,time) values, then keep WV fixed (ie do not fit)
settings.iFixO3_NoFit = -1;            %% -1 : do not fix O3 to ERA anomaly O3(z,time) values, then fit for O3
                                       %% +1 : do     fix O3 to ERA anomaly O3(z,time) values, then keep O3 fixed (ie do not fit)
                                       %% +0 : do     fix O3 to zero        O3(z,time) values, then keep O3 fixed (ie do not fit)

settings.iFixTG_NoFit = -1;            %% -1 means retrieve all trace gases [CO2 N2O CH4 CFC11 CFC12]
                                       %% eg [4 5] means do not do CFC11,CFC12
                                       %% eg [4] means do not do CFC11
settings.UMBCvsERAjac = -1;            %% do not adjust jacobian, based on handful of clear sky retrieval days
                                       %% [+1] means to adjust
settings.iLatX = 11;                   %% do latbins 01-11,54-64 using one set of params, 12-53 using another

topts_allowedparams = [{'ocb_set'},{'numchan'},{'chan_LW_SW'},{'iChSet'},{'set_tracegas'},{'offsetrates'},{'set_era5_cmip6_airsL3'},{'set_era5_cmip6_airsL3_WV_T_O3'},...
        	       {'addco2jacs'},{'obs_corr_matrix'},{'invtype'},{'tie_sst_lowestlayer'},{'iNlays_retrieve'},...
                       {'descORasc'},{'dataset'},{'iXJac'},{'co2lays'},{'iDoStrowFiniteJac'},...
                       {'iFixTz_NoFit'},{'iFixWV_NoFit'},{'iFixO3_NoFit'},{'iFixTG_NoFit'},{'iaSequential'}...
                       {'resetnorm2one'},{'UMBCvsERAjac'},{'iSergioCO2'},{'iAdjLowerAtmWVfrac'},{'iVersQRenorm'},{'rCoupleT_WV'},{'iLatX'}];

%disp('settings before')
%settings

if narginS == 3
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
   if (length(intersect(topts_allowedparams,optvar{i})) == 1)
     eval(sprintf('settings.%s = topts.%s;', optvar{i}, optvar{i}));
     str = ['junk = topts.' optvar{i} ';']; eval(str)
     fprintf(1,'  check_settings.m : will use user setting %30s = %4i \n',optvar{i},junk(1));
   else
     fprintf(1,'topts param not in allowed list ... %s \n',optvar{i});
     error('quitting ');
   end
 end
end

%disp('settings after')
%settings
