load convert_sergio_clearskygrid_obsonly_Q16.mat

plot(unique(X))
plot(unique(Y))

addpath /home/sergio/MATLABCODE/COLORMAP
figure(1); pcolor(unique(X),unique(Y),squeeze(b_asc(:,:,1520))'); colorbar; colormap(usa2); caxis([-0.25 +0.25]); shading flat
figure(2); pcolor(unique(X),unique(Y),squeeze(b_desc(:,:,1520))'); colorbar; colormap(usa2); caxis([-0.25 +0.25]); shading flat

addpath /home/sergio/MATLABCODE/PLOTTER
XX = X'; YY = Y'; figure(1);  z = squeeze(b_asc(:,:,1520))'; z = z(:); scatter_coast(XX(:),YY(:),50,z); title('dBT1231/dt ASC');
XX = X'; YY = Y'; figure(2);  z = squeeze(b_desc(:,:,1520))'; z = z(:); scatter_coast(XX(:),YY(:),50,z); title('dBT1231/dt DESC');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); scatter_coast(X(:),Y(:),50,reshape(squeeze(airs_noiseTtrue(:,:,1520)),1,4608)/sqrt(120)); title('AIRS NeDTtrue /sqrt(120)')
figure(4); scatter_coast(X(:),Y(:),50,reshape(squeeze(b_err_asc(:,:,1520)),1,4608)); title('unc b\_err\_asc')
figure(5); scatter_coast(X(:),Y(:),50,reshape(squeeze(b_err_desc(:,:,1520)),1,4608)); title('unc b\_err\_desc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see eg ../AIRS_new_clear_scan_August2019_AMT2020PAPER/get_rates.m
nedt_40_clear_latbins = load('../AIRS_new_clear_scan_August2019_AMT2020PAPER/btn_avg.mat');
figure(6); plot(nanmean(nedt_40_clear_latbins.btn_avg(:,:,1520),2)); title('Latbin40clr Nedt 1231 cm-1')

figure(7);
addpath /home/sergio/MATLABCODE/PLOTTER
latbin40 = equal_area_spherical_bands(20);
latbinx40 = 0.5*(latbin40(1:end-1) + latbin40(2:end));
plot(latbinx40,nanmean(nedt_40_clear_latbins.btn_avg(:,:,1520),2),'b.-',unique(YY),mean(b_err_desc(:,:,1520),1),'gx-',...
     unique(YY),mean(airs_noiseTtrue(:,:,1520),1)/sqrt(120),'ro-');
hl = legend('Latbin40 AMT paper','b\_err\_desc AIRS STM','true AIRS neDT/{\surd}(120)','location','best','fontsize',8,'Interpreter','tex');

title('{\itAe}^{-\alpha\itt}sin\beta{\itt} {\surd} \alpha<<\beta')
