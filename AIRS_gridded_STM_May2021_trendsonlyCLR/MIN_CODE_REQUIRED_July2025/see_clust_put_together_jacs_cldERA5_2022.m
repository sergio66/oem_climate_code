%% see clust_put_together_jacs_clrERA5.m
JOB = iLatBin;

%%%%%%%% ORIG CODE %%%%%%%
%% recall log10(X) = log(X)/log(10)
%% kcarta does Q d(BT)/dQ == dBT/d(logQ) = dBT/dQ/Q = Q d(BT)/dQ
%% but if we change to log10 then dBT/dlog10(Q) = dBT/d(log(Q)/log(10)) = log10 dBT/dlog(Q) = log10 Q dBT/dQ = log10 dBT/dlog(Q)
%% disp('WARNING here we use log10 jacs ie multiply gas jacs by log_e_(10) ie jac --> jac * log(10) = 2.3026')
%% disp('WARNING here we use log10 jacs ie multiply gas jacs by log_e_(10) ie jac --> jac * log(10) = 2.3026')
%% disp('WARNING here we use log10 jacs ie multiply gas jacs by log_e_(10) ie jac --> jac * log(10) = 2.3026')
%%%%%%%% ORIG CODE %%%%%%%

miaow = load('sarta_chans_for_l1c.mat');
ind2834to2645 = miaow.ichan;

if ~exist('iOldORNew')
  iOldORNew = +9;  %% the 19 year ERA5 cld, Q09
  iOldORNew = +5;  %% the 20 year ERA5 clr (Q0.97-->1) Q05
end

if iOldORNew == 9
  %% see ~sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Nov2022_startSept2002_endAug2022_trendsonly_cldy_Q09/
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_cld_subjacLatBin_kCARTA_ERA5_20yr_CLD_Q09_'           num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(JOB,'%02i') '.mat'];

  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_cld_subjacLatBin_kCARTA_ERA5_20yr_CLD_Q09_v2_'           num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_v2_' num2str(JOB,'%02i') '.mat'];

elseif iOldORNew == 5
  %% see ~sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Nov2022_startSept2002_endAug2022_trendsonly_cldy_Q09/
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_kCARTA_ERA5_20yr_'           num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_' num2str(JOB,'%02i') '.mat'];
else
  error('unknown iOldORNew')
end

%sarta = load(SARTAjac);

if iOldORNew == 9
  [h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/2012/FixedNAN/all4608_era5_full12months_Qcumulative09.rtp');
  thedir0 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Nov2022_startSept2002_endAug2022_trendsonly_cldy_Q09//';
  fprintf(1,'get_jac_fast --> see_clust_put_together_jacs_clrERA5_2022.m --> iOldORNew == 2022 (ERA5) JOB = %2i   \n',JOB);
elseif iOldORNew == 5
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp');
  thedir0 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Mar2023_startSept2002_endAug2022_trendsonly/';
  fprintf(1,'get_jac_fast --> see_clust_put_together_jacs_clrERA5_2022.m --> iOldORNew == 2022 (ERA5) JOB = %2i   \n',JOB);
end

lps = compute_lapse_rate(h,p);
profilejunk.nlays = p.nlevs(:,(iLatBin-1)*72+iLonBin)-1;
profilejunk.plays = plevs2plays(p.plevs(:,(iLatBin-1)*72+iLonBin));
profilejunk.ptemp = p.ptemp(:,(iLatBin-1)*72+iLonBin);
profilejunk.gas_1 = p.gas_1(:,(iLatBin-1)*72+iLonBin);
profilejunk.gas_3 = p.gas_3(:,(iLatBin-1)*72+iLonBin);
profilejunk.stemp = p.stemp((iLatBin-1)*72+iLonBin);
profilejunk.spres = p.spres((iLatBin-1)*72+iLonBin);
profilejunk.lps_tropoapauseP   = lps.trp_pHI((iLatBin-1)*72+iLonBin);
profilejunk.lps_tropoapauseind = lps.trp_ind((iLatBin-1)*72+iLonBin);

%% do the iNlaysavg later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor_log10 = log(10); %% this changes gas jac scaling to log10
factor_log10 = 1.0;     %% this keeps   gas jac scaling to loge

ind_lat_junk = JOB;
ind_subset_junk = (1:72) + (ind_lat_junk-1)*72;

for lon = iLonBin

  ind_subset_junk_ii = ind_subset_junk(lon);
  if iOldORNew == 9
    frad0 = [thedir0 '/AllDemJacsCld_Q09/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '.mat'];
    fz    = [thedir0 '/AllDemJacsCld_Q09/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_jac.mat'];
    fcol  = [thedir0 '/AllDemJacsCld_Q09/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_coljac.mat'];
  elseif iOldORNew == 5
    frad0 = [thedir0 '/AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '.mat'];
    fz    = [thedir0 '/AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_jac.mat'];
    fcol  = [thedir0 '/AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_coljac.mat'];
  end

  fprintf(1,'in see_clust_put_together_jacs_cldERA5_2022.m : iOldORNew  = %2i \n',iOldORNew);
  fprintf(1,'  frad0 = %s \n',frad0);
  fprintf(1,'  fz    = %s \n',fz);
  fprintf(1,'  fcol  = %s \n',fcol);

  iaIndices = ind_subset_junk_ii;

  arad0 = load(frad0);
  az    = load(fz);
  acol  = load(fcol); acol.rKc(:,1:6) = acol.rKc(:,1:6) * factor_log10;
  
  trad0 = rad2bt(az.fKc,arad0.rKc);
  %tcol  = rad2bt(az.fKc,acol.rKc);

  clear aout
  aout.fKc = az.fKc;

  %% kcarta nml = 2 4 5 6 51 52 T ST
  aout.jac(:,1) = acol.rKc(:,1);     %%% CO2
  colo3 = acol.rKc(:,2);             %%% (O3 is ind2, can compare againt the sum(o3jac(z))
  aout.jac(:,2) = acol.rKc(:,3);     %%% N2O (O3 is ind2)
  aout.jac(:,3) = acol.rKc(:,4);     %%% CH4
  aout.jac(:,4) = acol.rKc(:,5);     %%% CFC11
  aout.jac(:,5) = acol.rKc(:,6);     %%% CFC12
  aout.jac(:,6) = acol.rKc(:,8);     %%% ST
  
  [~,numlays] = size(az.rKc);
  numlays = (numlays-4)/4;
  
  figure(1); plot(aout.fKc,aout.jac); 
  
  ind = (1:numlays);
  
  iPlot = -1;

  ix = ind + numlays*0; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*factor_log10;  %% WV
    [an(1) an(end) ix(1) ix(end)]; 
    wvind = fliplr(an);
  if iPlot > 0
    figure(2); plot(aout.fKc,sum(az.rKc(:,ix),2)*factor_log10,'r',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacWV(1:numlays,:,lon)),1),'g'); title('WV')
     hl = legend('sum(Kc(1:97))','sarta','location','best','fontsize',10);
  end

  ix = ind + numlays*2; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix);  %% T
    [an(1) an(end) ix(1) ix(end)];
    tzind = fliplr(an);
  if iPlot > 0
    figure(3); plot(aout.fKc,sum(az.rKc(:,ix),2),'r',aout.fKc,acol.rKc(:,7),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacT(1:numlays,:,lon)),1),'g'); title('T')
%    figure(3); plot(aout.fKc,sum(az.rKc(:,ix),2),'r',aout.fKc,(tcol(:,7)-trad0),'b.-',...
%                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacT(1:numlays,:,lon)),1),'g'); title('T')
     hl = legend('sum(Kc(1:97))','col Kc','sarta','location','best','fontsize',10);
  end

   %% except we do not do col O3 in COL CLR jacs
  ix = ind + numlays*1; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*factor_log10;  %% O3
    [an(1) an(end) ix(1) ix(end)];
    o3ind = fliplr(an);
%{
    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*factor_log10,'r',aout.fKc,acol.rKc(:,2),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacO3(1:numlays,:,lon)),1),'g'); title('O3')
%    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*factor_log10,'r',aout.fKc,(tcol(:,2)-trad0)*factor_log10,'b.-',...
%                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacO3(1:numlays,:,lon)),1),'g'); title('O3')
     hl = legend('sum(Kc(1:97))','col Kc','sarta','location','best','fontsize',10);
%}

  ixx = numlays*4 + 1;  %% surface temp
  if iPlot > 0
    figure(1); plot(aout.fKc,aout.jac,aout.fKc,az.rKc(:,ixx),'b');   
    figure(5); plot(aout.fKc,az.rKc(:,ixx),'ro-',aout.fKc,aout.jac(:,6),'bx-',...
                  aout.fKc(ind2834to2645),squeeze(sarta.subjac.jacST(:,lon)),'g.');  title('SurfT')
     hl = legend('Kc','for retr','sarta','location','best','fontsize',10);
    figure(6); plot(aout.fKc,aout.jac(:,1),'b',...
                  aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacCO2z(1:numlays,:,lon)),1),'g'); title('CO2')
     hl = legend('for retr','sarta','location','best','fontsize',10);
  end

  aout.fKc = aout.fKc(ind2834to2645);
  aout.jac = aout.jac(ind2834to2645,:);
  colo3    = colo3(ind2834to2645);

  kcarta.subjac.coljacCO2   = aout.jac(:,1);
  kcarta.subjac.coljacN2O   = aout.jac(:,2);
  %kcarta.subjac.coljacCO   = aout.jac(:,X);
  kcarta.subjac.coljacCH4   = aout.jac(:,3);
  kcarta.subjac.coljacCFC11 = aout.jac(:,4);
  kcarta.subjac.coljacCFC12 = aout.jac(:,5);
  %  kcarta.subjac.jacCld1   = sarta.subjac.jacCld1;
  %  kcarta.subjac.jacCld2   = sarta.subjac.jacCld2;
  %  kcarta.subjac.jacSze1   = sarta.subjac.jacSze1;
  %  kcarta.subjac.jacSze2   = sarta.subjac.jacSze2;
  %  kcarta.subjac.jacTop1   = sarta.subjac.jacTop1;
  %  kcarta.subjac.jacTop2   = sarta.subjac.jacTop2;
  kcarta.subjac.jacST     = aout.jac(:,6);
  kcarta.subjac.jacT(1:numlays,:)  = aout.jac(:,tzind)';
  kcarta.subjac.jacWV(1:numlays,:) = aout.jac(:,wvind)';
  kcarta.subjac.jacO3(1:numlays,:) = aout.jac(:,o3ind)';
  kcarta.subjac.jacCO2z(1,:)       = aout.jac(:,1);

end
fprintf(1,'\n')

%kcarta.subjac.rlon = sarta.subjac.rlon(iLonBin);
%kcarta.subjac.rlat = sarta.subjac.rlat(iLonBin);
%kcarta.subjac.stemp = sarta.subjac.stemp(iLonBin);
%kcarta.subjac.nlevs = sarta.subjac.nlevs(iLonBin);
%kcarta.subjac.ppmv2 = sarta.subjac.ppmv2(iLonBin);
%kcarta.subjac.ppmv4 = sarta.subjac.ppmv4(iLonBin);
%kcarta.subjac.ppmv6 = sarta.subjac.ppmv6(iLonBin);

[~,kcarta.subjac.ppmv2] = layers2ppmv(h,p,iWhichLatLonBin,2);
[~,kcarta.subjac.ppmv4] = layers2ppmv(h,p,iWhichLatLonBin,4);
[~,kcarta.subjac.ppmv6] = layers2ppmv(h,p,iWhichLatLonBin,6);
kcarta.subjac.indices = iaIndices;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_ts_jac_fast = [];
%m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacCO2z/factor_log10];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCO2/factor_log10];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacN2O/factor_log10];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCH4/factor_log10];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCFC11/factor_log10];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCFC12/factor_log10];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacST];

m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacWV'];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacT'];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacO3'];

nlays = p.nlevs(iWhichLatLonBin)-1;
[~,nlaysx] = size(m_ts_jac_fast);
nlaysx = (nlaysx-6)/3;
if nlaysx ~= nlays
  fprintf(1,'<<<<<<<<<<<<< warning see_clust_put_together_jacs_clrERA5_2021.m nlaysx ~= nlays %3i %3i >>>>>>>>>>>>>> \n',nlaysx,nlays);
  nlays = nlaysx;
end

freq2645 = h.vchan;
