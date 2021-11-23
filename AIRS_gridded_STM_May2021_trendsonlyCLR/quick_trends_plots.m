addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/science/

if ~exist('b_obs') & ~exist('b_asc') & ~exist('b_desc')
  iInput = input('(+1) Allsky or (-1) ClrSky : ');
  if iInput == +1
    load xconvert_strowrates2oemrates_allskygrid_obsonly.mat
  else
    load convert_sergio_clearskygrid_obsonly_Q16.mat
    b_obs = b_desc;
  end  
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

iChan = input('Enter wavenumber to display (645:2780), -1 to stop) : ');
while iChan > 0
  iChan = find(h.vchan >= iChan,1);
  fprintf('showing h2645.iChan = %4i h2645.vchan = %9.4f \n',h.ichan(iChan),h.vchan(iChan))
  if iInput == +1
    figure(1); junk = squeeze(b_obs(:,iChan,1));
    scatter_coast(X(:),Y(:),50,junk); colorbar; colormap(usa2);
  else
    figure(2); junk = squeeze(b_obs(:,:,iChan));
    scatter_coast(X(:),Y(:),50,junk(:)); colorbar; colormap(usa2);
    cx = caxis; caxis([-0.1 +0.1]); colorbar;
  end
  iChan = input('Enter wavenumber to display (6045:2780), -1 to stop) : ');
end

