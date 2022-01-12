%thedir0 = '/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv2/LatBin' num2str(JOB,'%02d') '/';
%iLon = 1;
%thefilein  = [thedir0 '/' LonBin' num2str(iLon,'%02d') '/cfbins_LatBin' num2str(JOB,'%02d') '_LonBin' num2str(iLon,'%02d') '_V2.mat'];

addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/science/

% if ~exist('b_obs')
%   %load convert_strowrates2oemrates_allskygrid_obsonly.mat
%   load convert_sergio_clearskygrid_obsonly_Q16.mat
%   b_obs = b_desc;
% end

if ~exist('h')
  load h2645structure.mat
end

airs_noise = instr_chans2645('airs',2);

load latB64.mat
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
figure(1); scatter_coast(X(:),Y(:),50,landfrac); colorbar; title('landfrac');  caxis([0 1])

quants =  [0.0100 0.0200 0.0300 0.0400 0.0500 0.1000 0.2500 0.5000 0.7500 0.9000 0.9500 0.9600 0.9700 0.9800 0.9900 1]

b_asc = nan(72,64,2645);
b_desc = nan(72,64,2645);

iType = +01; %% strows Q16 18 year trends       Strow ran for me in Feb 2021, 2002/09 to 2020/08
iType = -01; %% sergio Q16 18 year TEST trends  I     ran for me in Aug 2021, 2002/09 to 2020/08, to check I get same results as Strow +01
iType = +02; %% sergio Q16 19 year trends       I     ran for me in Aug 2021, 2002/09 to 2021/07
iType = +03; %% sergio Extreme 19 year trends   I     ran for me in Aug 2021, 2002/09 to 2020/08
iType = -03; %% sergio Mean 19 year trends      I     ran for me in Aug 2021, 2002/09 to 2020/08
iType = +04; %% sergio Q16 19 year trends       I     ran for me in Aug 2021, 2002/09 to 2021/08 FULL

disp('Choices DataSet to use ')
disp('                       (+1) Strow  Quantile Mar 2021 2002/09 to 2020/08 Full 18 years');
disp('                       (-1) Sergio Quantile Aug 2021 2002/09 to 2020/08 Full 18 years');
disp('                        (2) Sergio Quantile Aug 2021 2002/09 to 2021/07 ');
disp('                        (4) Sergio Quantile Aug 2021 2002/09 to 2021/08 Full 19 years **** ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                        (3) Sergio Extreme Aug 2021 2002/09 to 2021/07 ');
disp('                       (-3) Sergio Mean Aug 2021 2002/09 to 2021/07 ');
disp(' <---------------------------------------------------------------------------------------> ')
iType = input('Enter DataSet to use (+1,-1,+2,+4   or +3,-3) : ');

if iType ~= 3
  iQuantile = 16;  %% hottest, used for AIRS STM May 221
  iQuantile = 08;
  iQuantile = 04;
  iQuantile = input('Enter iQuantile to make (1-16, 0 = avg) : ');
end

iKeepPlotting = -1;

for iLat = 1 : 64
  if iType == +1
    thedir0 = ['/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv1_fits/LatBin' num2str(iLat,'%02d') '/'];   %% STROW STUFF March 2021
  elseif iType == -1
    % see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 2
    % see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 3
    % see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == -3
    % see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon_v3/Mean/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 4
    % see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
  end

  for iLon = 1 : 72
    fprintf(1,'+')
    if iType == +1
      %% full 18 year by Strow
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1.mat'];
    elseif iType == -1
      %% full 18 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps412.mat'];
    elseif iType == 2
      %% paritally incomplete 19 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps429.mat'];
    elseif iType == 3
      %% paritally incomplete 19 year by Sergio
      thefilein = [thedir0 'LonBin' num2str(iLon,'%02d') '/extreme_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps429.mat'];
    elseif iType == -3
      %% paritally incomplete 19 year by Sergio
      thefilein = [thedir0 'LonBin' num2str(iLon,'%02d') '/mean_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps429.mat'];
    elseif iType == 4
      %% full 19 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps433.mat'];
    end

    iBoo = (iLat-1)*72 + iLon;
    if ~exist(thefilein)
      fprintf(1,'oops iLat/iLon = %2i %2i tilenumber %4i file %s DNE \n',iLat,iLon,iBoo,thefilein);
    end
    x = load(thefilein);

    if abs(iType) ~= 3
      %% QUANTILES QUANTILES QUANTILES
      %% these are rads
      %b_asc(iLon,iLat,:) = x.b_asc(:,iQuantile,2);
      %b_desc(iLon,iLat,:) = x.b_desc(:,iQuantile,2);
      %b_err_asc(iLon,iLat,:) = x.berr_asc(:,iQuantile,2);
      %b_err_desc(iLon,iLat,:) = x.berr_desc(:,iQuantile,2);
  
      %% these are bt
      b_asc(iLon,iLat,:) = x.dbt_asc(:,iQuantile);
      b_desc(iLon,iLat,:) = x.dbt_desc(:,iQuantile);
      b_err_asc(iLon,iLat,:) = x.dbt_err_asc(:,iQuantile);
      b_err_desc(iLon,iLat,:) = x.dbt_err_desc(:,iQuantile);
      lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc(:,iQuantile);
      lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc(:,iQuantile);
  
      mean_rad(iLon,iLat,:) = squeeze(x.b_desc(:,iQuantile,1)); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
      mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,squeeze(x.b_desc(:,iQuantile,1)));
  
      junk = squeeze(mean_BT(iLon,iLat,:));
      junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      airs_noiseTtrue(iLon,iLat,:) = junk;

    elseif iType == +3
      %% these are extremes      
      b_asc(iLon,iLat,:) = x.dbt_asc;
      b_desc(iLon,iLat,:) = x.dbt_desc;
      b_err_asc(iLon,iLat,:) = x.dbt_err_asc;
      b_err_desc(iLon,iLat,:) = x.dbt_err_desc;
      lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc';
      lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc';
  
      mean_rad(iLon,iLat,:) = x.b_desc(:,1); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
      mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,mean_rad(iLon,iLat,:));
  
      junk = squeeze(mean_BT(iLon,iLat,:));
      junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      airs_noiseTtrue(iLon,iLat,:) = junk;

    elseif iType == -3
      %% these are means      
      b_asc(iLon,iLat,:) = x.dbt_asc;
      b_desc(iLon,iLat,:) = x.dbt_desc;
      b_err_asc(iLon,iLat,:) = x.dbt_err_asc;
      b_err_desc(iLon,iLat,:) = x.dbt_err_desc;
      lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc';
      lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc';
  
      mean_rad(iLon,iLat,:) = x.b_desc(:,1); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
      mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,mean_rad(iLon,iLat,:));
  
      junk = squeeze(mean_BT(iLon,iLat,:));
      junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      airs_noiseTtrue(iLon,iLat,:) = junk;

    end

  end
  fprintf(1,'%3i out of 64 \n',iLat);

  if iKeepPlotting > 0
    junk = squeeze(b_desc(:,:,1520));
    scatter_coast(X(:),Y(:),50,junk(:)); pause(0.1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iType == 1
  saver = ['save convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
elseif iType == -1
  saver = ['save iType_-1_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
elseif iType == 2
  saver = ['save iType_2_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
elseif iType == 4
  saver = ['save iType_4_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
%%
%%
elseif iType == 3
  saver = ['save iType_3_extreme_convert_sergio_clearskygrid_obsonly.mat b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
elseif iType == -3
  saver = ['save iType_-3_mean_convert_sergio_clearskygrid_obsonly.mat b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
end
eval(saver);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
[g2378,g2645] = compare_goodchans_2378_2645();   %% shows g is index into f-ABCD and NOT chanID

figure(1)
junk = squeeze(b_desc(:,:,1520));
scatter_coast(X(:),Y(:),50,junk(:)); pause(0.1);
caxis([-0.2 +0.2]);; colormap(usa2)
title(['Quantile ' num2str(iQuantile,'%02d')])

figure(2); clf
junk = rad2bt(1231,squeeze(mean_rad(:,:,1520)));
scatter_coast(X(:),Y(:),50,junk(:));
colormap jet
title(['Quantile ' num2str(iQuantile,'%02d')])

figure(2); clf
junk = squeeze(mean_BT(:,:,1520));
scatter_coast(X(:),Y(:),50,junk(:));
colormap jet
title(['Quantile ' num2str(iQuantile,'%02d')])

figure(3); clf
junk = squeeze(airs_noiseTtrue(:,:,1520));
junk = reshape(airs_noiseTtrue,72*64,2645)'/sqrt(120);    %% need sqrt(N) from about 12000 obs/tile/16 days .. so 1% of this is 120 ... noise goes down by sqrt(N)
plot(h.vchan,nanmean(reshape(b_err_desc,72*64,2645)',2),h.vchan(g2645),nanmean(junk(g2645,:),2));  
  ylim([0 0.3]); 
  hl = legend('from b_err','from 1/sqrt(N)','location','best','fontsize',10);
title(['Quantile ' num2str(iQuantile,'%02d')])

figure(4); clf
junkA = reshape(b_asc,4608,2645);
junkD = reshape(b_desc,4608,2645);
plot(h.vchan,nanmean(junkA,1),h.vchan,nanmean(junkD,1),h.vchan,nanstd(junkA,1),'--',h.vchan,nanstd(junkD,1),'--')
  ylim([-0.1 +0.1]*0.75)
  xlim([640 1640])
  hl = legend('mean ratesA','mean ratesD','std ratesA','std ratesD','location','best','fontsize',10);
title(['Quantile ' num2str(iQuantile,'%02d')])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath /home/sergio/KCARTA/MATLAB
plot(h.vchan,reshape(b_err_desc,72*64,2645)')
plot(h.vchan,reshape(mean_BT,72*64,2645)')
airs_noise = instr_chans2645('airs',2);
for iLat = 1 : 64
  for iLon = 1 : 72
    junk = squeeze(mean_BT(iLon,iLat,:));
    junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
    airs_noiseTtrue(iLon,iLat,:) = junk;
  end
end

addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
[g2378,g2645] = compare_goodchans_2378_2645();   %% shows g is index into f-ABCD and NOT chanID
junk = reshape(airs_noiseTtrue,72*64,2645)'/sqrt(120);   %% about 12000 obs/tile/16 days .. so 1% of this is 120 ... noise goes down by sqrt(N)
plot(h.vchan(g2645),junk(g2645,:))
ylim([0 0.3])

plot(h.vchan,nanmean(reshape(b_err_desc,72*64,2645)',2),h.vchan(g2645),nanmean(junk(g2645,:),2)); 
  ylim([0 0.3]); 
  hl = legend('from b_err','from 1/sqrt(N)','location','best','fontsize',10);
%}
