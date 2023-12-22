addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

load llsmap5

if ~exist('b_desc')
  load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat
end

do_XX_YY_from_X_Y

figure(1); pcolor(squeeze(b_desc(:,:,1520))'); caxis([-1 +1]*0.15); colormap(usa2); title('BT 1231 trend : DESC'); shading interp
figure(2); pcolor(squeeze(b_asc(:,:,1520))'); caxis([-1 +1]*0.15); colormap(usa2); title('BT 1231 trend : ASC'); shading interp
figure(3); pcolor(0.5*squeeze(b_asc(:,:,1520)+b_desc(:,:,1520))'); caxis([-1 +1]*0.15); colormap(usa2); title('BT 1231 trend : 0.5*(DESC+ASC)'); shading interp

figure(1); scatter_coast(X,Y,100,squeeze(b_desc(:,:,1520))); caxis([-1 +1]*0.15); colormap(llsmap5); title('BT 1231 trend : DESC'); shading interp
figure(2); scatter_coast(X,Y,100,squeeze(b_asc(:,:,1520))); caxis([-1 +1]*0.15); colormap(llsmap5); title('BT 1231 trend : ASC'); shading interp
figure(3); scatter_coast(X,Y,100,0.5*squeeze(b_asc(:,:,1520)+b_desc(:,:,1520))); caxis([-1 +1]*0.15); colormap(llsmap5); title('BT 1231 trend : 0.5*(DESC+ASC)'); shading interp
