%% see clust_put_together_jacs_clr.m
JOB = iLatBin;

miaow = load('sarta_chans_for_l1c.mat');
ind2834to2645 = miaow.ichan;

iOldORNew = +1;
if iOldORNew < 0
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  %%% NEED TO REDO
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin' num2str(JOB,'%02i') '.mat'];
else
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_newSARTA_' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_newSARTA_' num2str(JOB,'%02i') '.mat'];
end

%sarta = load(SARTAjac);

%[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019_palts.rtp');  %% already has palts

%[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp');  %% already has palts
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
%plot(1:72,p.stemp(sarta.subjac.indices),'b-',1:72,sarta.subjac.stemp,'r.-')
%pause(1); 

ind_lat_junk = JOB;
ind_subset_junk = (1:72) + (ind_lat_junk-1)*72;

thedir0 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Oct2020_startSept2002_trendsonly/';

for lon = iLonBin

  ind_subset_junk_ii = ind_subset_junk(lon);
  frad0 = [thedir0 '/AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '.mat'];
  fz    = [thedir0 '/AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_jac.mat'];
  fcol  = [thedir0 '/AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_coljac.mat'];
  iaIndices = ind_subset_junk_ii;

  arad0 = load(frad0);
  az    = load(fz);
  acol  = load(fcol); acol.rKc(:,1:6) = acol.rKc(:,1:6) * log(10);
  
  trad0 = rad2bt(az.fKc,arad0.rKc);
  %tcol  = rad2bt(az.fKc,acol.rKc);

  clear aout
  aout.fKc = az.fKc;

  %% kcarta nml = 2 4 5 6 51 52 T ST
  aout.jac(:,1) = acol.rKc(:,1);     %%% CO2
  aout.jac(:,2) = acol.rKc(:,2);     %%% N2O
  aout.jac(:,3) = acol.rKc(:,4);     %%% CH4
  aout.jac(:,4) = acol.rKc(:,5);     %%% CFC11
  aout.jac(:,5) = acol.rKc(:,6);     %%% CFC12
  aout.jac(:,6) = acol.rKc(:,8);     %%% ST
  
  [~,numlays] = size(az.rKc);
  numlays = (numlays-4)/4;
  
  figure(1); plot(aout.fKc,aout.jac); 
  
  ind = (1:numlays);
  
  iPlot = -1;

  ix = ind + numlays*0; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*log(10);  %% WV
    [an(1) an(end) ix(1) ix(end)]; 
    wvind = fliplr(an);
  if iPlot > 0
    figure(2); plot(aout.fKc,sum(az.rKc(:,ix),2)*log(10),'r',...
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
  ix = ind + numlays*1; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*log(10);  %% O3
    [an(1) an(end) ix(1) ix(end)];
    o3ind = fliplr(an);
%{
    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*log(10),'r',aout.fKc,acol.rKc(:,2),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacO3(1:numlays,:,lon)),1),'g'); title('O3')
%    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*log(10),'r',aout.fKc,(tcol(:,2)-trad0)*log(10),'b.-',...
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
%m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacCO2z/log(10)];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCO2/log(10)];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacN2O/log(10)];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCH4/log(10)];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCFC11/log(10)];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.coljacCFC12/log(10)];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacST];

nlays = p.nlevs(iWhichLatLonBin)-1;
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacWV'];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacT'];
m_ts_jac_fast = [m_ts_jac_fast kcarta.subjac.jacO3'];

freq2645 = h.vchan;
