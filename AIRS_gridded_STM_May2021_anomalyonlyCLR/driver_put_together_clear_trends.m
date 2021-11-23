%thedir0 = '/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv2/LatBin' num2str(JOB,'%02d') '/';
%iLon = 1;
%thefilein  = [thedir0 '/' LonBin' num2str(iLon,'%02d') '/cfbins_LatBin' num2str(JOB,'%02d') '_LonBin' num2str(iLon,'%02d') '_V2.mat'];

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/science/

if ~exist('b_obs')
  %load convert_strowrates2oemrates_allskygrid_obsonly.mat
  load convert_sergio_clearskygrid_obsonly.mat
  b_obs = b_desc;
end

if ~exist('h')
  load h2645structure.mat
end

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

for iLat = 1 : 64
  thedir0 = ['/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv1_fits/LatBin' num2str(iLat,'%02d') '/'];
  for iLon = 1 : 72
    fprintf(1,'+')
    thefilein  = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1.mat'];
    x = load(thefilein);

    %% these are rads
    %b_asc(iLon,iLat,:) = x.b_asc(:,16,2);
    %b_desc(iLon,iLat,:) = x.b_desc(:,16,2);
    %b_err_asc(iLon,iLat,:) = x.berr_asc(:,16,2);
    %b_err_desc(iLon,iLat,:) = x.berr_desc(:,16,2);

    %% these are bt
    b_asc(iLon,iLat,:) = x.dbt_asc(:,16);
    b_desc(iLon,iLat,:) = x.dbt_desc(:,16);
    b_err_asc(iLon,iLat,:) = x.dbt_err_asc(:,16);
    b_err_desc(iLon,iLat,:) = x.dbt_err_desc(:,16);
    lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc(:,16);
    lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc(:,16);
  end
  fprintf(1,'%3i out of 64 \n',iLat);

  junk = squeeze(b_desc(:,:,1520));
  scatter_coast(X(:),Y(:),50,junk(:)); pause(0.1);
end
junk = squeeze(b_desc(:,:,1520));
scatter_coast(X(:),Y(:),50,junk(:)); pause(0.1);
caxis([-0.2 +0.2]);

save convert_sergio_clearskygrid_obsonly.mat b_* X Y landfrac salti h lagcor*
