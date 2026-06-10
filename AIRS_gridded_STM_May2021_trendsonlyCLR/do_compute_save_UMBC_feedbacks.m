feedbacknameUMBC = 'ThisDoesNotExist';
if ~exist('umbc_spectral_olr') 
  if exist(feedbacknameUMBC)
    iUMBCexist = +1;
    fprintf(1,'for %2i years the UMBC OLR ecRad file %s exists \n',iNumYears,feedbacknameUMBC)
    lser = ['!ls -lt ' feedbacknameUMBC];
    eval(lser)
    junk = input('re-do the UMBC OLR feedbacks??? (-1 [default] /+1) : ');    
    if junk > 0
      umbc_spectral_olr = struct;    %% so it has no fields
      umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,-1,rlat65,rlon73,'UMBC');
      %%% compute_feedbacks_umbc_ecRad   ; pause(0.1)
    else
      fprintf(1,'  reading in %s \n',feedbacknameUMBC)
      junk = load(feedbacknameUMBC,'umbc_spectral_olr','results');
      umbc_spectral_olr = junk.umbc_spectral_olr;
      umbc_spectral_olr.stemptrend = junk.results(:,6)';
      clear junk
    end
  else
    iTest100profiles = -1;  %% default
    if iTest100profiles > 0
      do_ecrad_UMBC_Test100profiles
    end
    
    disp('doing compare_two_structures(p,pMean17years) since we send "p" in ')
    compare_two_structures(p,pMean17years);
    if h.pfields == 1 | ~isfield(p,'rcalc')
      disp('making p.rcalc ....')
      pcopy = p;
      hcopy = h;      
      klayers = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_80km/klayersV205_0_80km/BinV221/klayers_airs';
      sarta   = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_apr26_H2024';
      rtpwrite('olrjunk.op.rtp',h,ha,p,pa);
      sartaer = ['!time ' sarta ' fin=olrjunk.op.rtp fout=olrjunk.rp.rtp'];
      eval(sartaer)
      [h,~,p,~] = rtpread('olrjunk.rp.rtp');
      rmer = ['!/bin/rm olrjunk.op.rtp olrjunk.rp.rtp'];
      eval(rmer)
    end
    umbc_spectral_olr       = struct;    %% so it has no fields    
    umbc_spectral_olr       = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,-1,rlat65,rlon73,-1,'UMBC');
    umbc_spectral_delta_olr = struct;    %% so it has no fields        
    umbc_spectral_delta_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_delta_olr,8888,rlat65,rlon73,-1,'UMBC');    
    %%% compute_feedbacks_umbc_ecRad   ; pause(0.1)
  end
end

junk = input('save the UMBC OLR feedbacks??? (-1 [default] /+1) : ');
if length(junk) == 0
  junk = -1;
end
if junk > 0
  %{
  %%%% this is how I made  olr_feedbacks_globalSST_AIRSL3_ERA5_CMIP6_save.mat olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat %%%%
  if iSwap_ERA_2012_08_15 < 0
    feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];
  else
    feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_swap_profile.mat'];
  end
  if iRaw_or_Unc == -1
    feedbacknameUMBC = feedbacknameUMBC(1:end-4);
    feedbacknameUMBC = [feedbacknameUMBC '_unc_factor' num2str(maxratio,'%0.2f') '.mat'];
  end
  saver = ['save ' feedbackname ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays   airsL3_spectral_olr era5_spectral_olr cmip6_spectral_olr airsL3 era5 cmip6 strUMBC iNumYears'];
  saver = ['save ' feedbackname ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays   airsL3_spectral_olr era5_spectral_olr cmip6_spectral_olr airsL3 era5 cmip6 strUMBC iNumYears'];
  %%%% this is how I made  olr_feedbacks_globalSST_AIRSL3_ERA5_CMIP6_save.mat olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat %%%%
  %}
  
  if junk > 0
    junk2 = +1;
    if exist(feedbacknameUMBC)
      lser = ['!ls -lth ' feedbacknameUMBC];
      eval(lser);
      junk2 = input('file already exists, overwrite (-1 default/+1) : ');
      if length(junk2) == 0
        junk2 = -1;
      end
    end
    if ~exist('strUMBC')
      strUMBC = 'tempStrUMBC';
    end    
    saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr umbc_spectral_delta_olr results resultsWV resultsT resultsO3 pavg plays strUMBC iNumYears'];  %% if you only want to save UMBC
    if junk2 > 0
      fprintf(1,'saving to %s \n',feedbacknameUMBC);
      eval(saver);
      %%%% mver = ['!mv ThisDoesNotExist.mat umbc_only_OLR_feedbacks_23years_try1.mat']; eveal(mver)
    else
      fprintf(1,'this already exists %s not saving \n',feedbacknameUMBC);
    end
  end
end
