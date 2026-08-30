disp('see the AK test vs Bill Irion in /home/sergio/MATLABCODE_Git/AveragingKernel_DOF/driver_compare_AK_JPL_vs_UMBC.m')

%load Output/Quantile03/test1000.mat
%load Output/Quantile03/test76.mat
load Output/Quantile03/test2000.mat

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN./playsD;
plays = flipud(plays);

iNumLay = 49;
iNumLay = length(jacobian.wvjaclays_used);
if length(jacobian.wvjaclays_used) == iNumLay
  for iii = 1 : length(jacobian.wvjaclays_used)
    %iavg = jacobian.wvjaclays_used{iNumLay}-6;
    iavg = jacobian.wvjaclays_used{iii}-jacobian.wvjaclays_offset;
    pavg(iii) = mean(plays(iavg));
  end
end

pmin = 0.1;
pmin = 1.0;
pmin = 10.0;
pmin = 80.0;
figure(1); clf; semilogy(oem.ak_ozone',pavg); ylim([pmin 1000]); title('OZ(p)'); set(gca,'ydir','reverse')
figure(2); clf; semilogy(oem.ak_water',pavg); ylim([pmin 1000]); title('WV(p)'); set(gca,'ydir','reverse')
figure(3); clf; semilogy(oem.ak_temp',pavg);  ylim([pmin 1000]); title('TZ(p)'); set(gca,'ydir','reverse')
figure(4); clf; semilogy(nanmean(oem.ak_ozone),pavg); ylim([pmin 1000]); title('OZ(p)'); set(gca,'ydir','reverse')
figure(5); clf; semilogy(nanmean(oem.ak_water),pavg); ylim([pmin 1000]); title('WV(p)'); set(gca,'ydir','reverse')
figure(6); clf; semilogy(nanmean(oem.ak_temp),pavg);  ylim([pmin 1000]); title('TZ(p)'); set(gca,'ydir','reverse')
figure(7); clf; semilogy(nanmean(oem.ak_ozone),pavg,'k',nanmean(oem.ak_water),pavg,'b',nanmean(oem.ak_temp),pavg,'r','linewidth',2); 
  ylim([pmin 1000]); title('ALL (p)'); set(gca,'ydir','reverse'); hl = legend('OZ','WV','T','location','best','fontsize',10); grid on;
figure(8); clf; semilogy(nansum(oem.ak_temp,2),pavg,'r',nansum(oem.ak_water,2),pavg,'b',nansum(oem.ak_ozone,2),pavg,'k','linewidth',2); ylim([pmin 1000]); 
  title('Row Sum of AK'); set(gca,'ydir','reverse'); legend('T','WV','O3','location','best','fontsize',10);
