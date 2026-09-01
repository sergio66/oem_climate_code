addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/science/

if ~exist('a16')
  a16 = load('convert_sergio_clearskygrid_obsonly_Q16.mat');
  a08 = load('convert_sergio_clearskygrid_obsonly_Q08.mat');
  a04 = load('convert_sergio_clearskygrid_obsonly_Q04.mat');
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

cx = 0.25;  cx = 0.10;
yoff = 0.5; yoff = 0.05;

iChan = input('Enter wavenumber to display (645:2780), -1 to stop) : ');
while iChan > 0
  iChan = find(h.vchan >= iChan,1);
  fprintf('showing h2645.iChan = %4i h2645.vchan = %9.4f \n',h.ichan(iChan),h.vchan(iChan))

  figure(1); junk = squeeze(a04.b_desc(:,:,iChan));
  scatter_coast(X(:),Y(:),50,junk(:)); colorbar; colormap(usa2); title(['Q04 ' num2str(h.vchan(iChan)) ' cm-1'])
  caxis([-cx +cx]); colorbar;

  figure(2); junk = squeeze(a08.b_desc(:,:,iChan));
  scatter_coast(X(:),Y(:),50,junk(:)); colorbar; colormap(usa2); title(['Q08 ' num2str(h.vchan(iChan)) ' cm-1'])
  caxis([-cx +cx]); colorbar;

  figure(3); junk = squeeze(a16.b_desc(:,:,iChan));
  scatter_coast(X(:),Y(:),50,junk(:)); colorbar; colormap(usa2); title(['Q16 ' num2str(h.vchan(iChan)) ' cm-1'])
  caxis([-cx +cx]); colorbar;

  btanom04 = squeeze(nanmean(a04.b_desc,1));
  btanom08 = squeeze(nanmean(a08.b_desc,1));
  btanom16 = squeeze(nanmean(a16.b_desc,1));
  tropics = find(abs(rlat) <= 30);
  figure(4);
  plot(h.vchan,nanmean(btanom04(tropics,:),1)-yoff,'b',h.vchan,nanmean(btanom08(tropics,:),1),'g',h.vchan,nanmean(btanom16(tropics,:),1)+yoff,'r','linewidth',2)
    xlim([645 1640]); grid; hl = legend('Q04','Q08','Q16','location','best','fontsize',10); 

  figure(5);
  plot(h.vchan,nanmean(btanom04(tropics,:),1)-yoff,'b',h.vchan,nanmean(btanom08(tropics,:),1),'g',h.vchan,nanmean(btanom16(tropics,:),1)+yoff,'r','linewidth',2)
    xlim([645 845]); grid; hl = legend('Q04','Q08','Q16','location','best','fontsize',10);

  figure(6);
  plot(h.vchan,nanmean(btanom04(tropics,:),1)-yoff,'b',h.vchan,nanmean(btanom08(tropics,:),1),'g',h.vchan,nanmean(btanom16(tropics,:),1)+yoff,'r','linewidth',2)
    xlim([1240 1640]); grid; hl = legend('Q04','Q08','Q16','location','best','fontsize',10);

  iChan = input('Enter wavenumber to display (6045:2780), -1 to stop) : ');
end

