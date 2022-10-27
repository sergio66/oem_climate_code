addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/TIME

if ~exist('woo')
  %woo = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Aug2014_12yr_desc.mat');
  woo = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Aug2021_19yr_desc.mat');
  [Y,X] = meshgrid(woo.save_lat64x72,woo.save_lon64x72);
  [salti, landfrac] = usgs_deg10_dem(Y,X);
  woo.landfrac = landfrac';
  woo.lat = Y';
  woo.lon = X';
end

ii = 0;
for yy = 2002 : 2021
  mmS = 1; mmE = 12;
  if yy == 2002
    mmS = 9;
  elseif yy == 2021
    mmE = 8;
  end
  for mm = mmS : mmE
    ii = ii + 1;
    yyx(ii) = yy;
    mmx(ii) = mm;
    ddx(ii) = 15;
  end
end

daysSince2002 = change2days(yyx,mmx,ddx,2002);

clear *fftspectrum*

%ii = 36; jj = 32;   fftspectrum_clr = fft(squeeze(woo.save64x72_clrolr(jj,ii,:)));
%plot(atan(imag(fftspectrum_clr)./real(fftspectrum_clr))*180/pi)
%semilogy(abs(fftspectrum_clr))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ylat = Y;

disp('Enter (-7) polar    L/O')
disp('      (-6) midlat   L/O')
disp('      (-5) tropical L/O')
disp('      (-4) polar land    (+4) polar ocean')
disp('      (-3) midlat land   (+3) midlat ocean')
disp('      (-2) tropical land (+2) tropical ocean')
disp('      (-1) land          (+1) ocean');
disp('      [0,default] ALL trends : ');
iAorOorL = input('Enter region : ');
if length(iAorOorL) == 0
  iAorOorL = 0;
end

clear maskLF
maskLF = nan(1,4608);
if iAorOorL == -7
  maskLF(abs(Ylat) > 60) = 1;
elseif iAorOorL == -6
  maskLF(abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == -5
  maskLF(abs(Ylat) < 30) = 1;
elseif iAorOorL == 0
  maskLF = ones(1,4608);
elseif iAorOorL == -1
  maskLF(landfrac == 1) = 1;
elseif iAorOorL == +1
  maskLF(landfrac == 0) = 1;
elseif iAorOorL == -2
  maskLF(landfrac == 1 & abs(Ylat) <= 30) = 1;
elseif iAorOorL == +2
  maskLF(landfrac == 0 & abs(Ylat) <= 30) = 1;
elseif iAorOorL == -3
  maskLF(landfrac == 1 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == +3
  maskLF(landfrac == 0 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == -4
  maskLF(landfrac == 1 & abs(Ylat) > 60) = 1;
elseif iAorOorL == +4
  maskLF(landfrac == 0 & abs(Ylat) > 60) = 1;
end
maskLFmatr = reshape(maskLF,72,64)';
mask = find(maskLF == 1);
figure(5); clf; pcolor(maskLFmatr); colorbar; title('maskLF')

[mm,nn,oo] = size(woo.save64x72_olr);
for ii = 1 : oo
  woo.xsave64x72_clrolr(:,:,ii) = squeeze(woo.save64x72_clrolr(:,:,ii)) .* maskLFmatr;  
  woo.xsave64x72_olr(:,:,ii)    = squeeze(woo.save64x72_olr(:,:,ii)) .* maskLFmatr;
  woo.xsave64x72_stemp(:,:,ii)  = squeeze(woo.save64x72_stemp(:,:,ii)) .* maskLFmatr;
  woo.xsave64x72_RHSurf(:,:,ii) = squeeze(woo.save64x72_RHSurf(:,:,ii)) .* maskLFmatr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
%% no detrending
for jj = 1 : 64
  wfftspectrum_clr(jj,:)    = fft(squeeze(nanmean(woo.xsave64x72_clrolr(jj,:,:),2)));
  wfftspectrum_cld(jj,:)    = fft(squeeze(nanmean(woo.xsave64x72_olr(jj,:,:),2)));
  wfftspectrum_stemp(jj,:)  = fft(squeeze(nanmean(woo.xsave64x72_stemp(jj,:,:),2)));
  wfftspectrum_RHSurf(jj,:) = fft(squeeze(nanmean(woo.xsave64x72_RHSurf(jj,:,:),2)));
end

iMax = length(wfftspectrum_cld);

%% fit for sin nw)
for jj = 1 : 64
  xfftspectrum_clr(jj,:)    = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_clrolr(jj,:,:),2)),1:iMax,1,4));
  xfftspectrum_cld(jj,:)    = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_olr(jj,:,:),2)),1:iMax,1,4));
  xfftspectrum_stemp(jj,:)  = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_stemp(jj,:,:),2)),1:iMax,1,4));
  xfftspectrum_RHSurf(jj,:) = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_RHSurf(jj,:,:),2)),1:iMax,1,4));
end

%% fit for sin 1/n w)
for jj = 1 : 64
  yfftspectrum_clr(jj,:)    = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_clrolr(jj,:,:),2)),1:iMax,-1,4));
  yfftspectrum_cld(jj,:)    = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_olr(jj,:,:),2)),1:iMax,-1,4));
  yfftspectrum_stemp(jj,:)  = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_stemp(jj,:,:),2)),1:iMax,-1,4));
  yfftspectrum_RHSurf(jj,:) = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_RHSurf(jj,:,:),2)),1:iMax,-1,4));
end

%% fit for sin nw + sin 1/n w
for jj = 1 : 64
  zfftspectrum_clr(jj,:)    = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_clrolr(jj,:,:),2)),1:iMax,0,4));
  zfftspectrum_cld(jj,:)    = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_olr(jj,:,:),2)),1:iMax,0,4));
  zfftspectrum_stemp(jj,:)  = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_stemp(jj,:,:),2)),1:iMax,0,4));
  zfftspectrum_RHSurf(jj,:) = fft(generic_compute_anomaly(daysSince2002,squeeze(nanmean(woo.xsave64x72_RHSurf(jj,:,:),2)),1:iMax,0,4));
end

iWhichFFT = input('Enter which fft (+10,default) = no detrend (1) = use sin nwt (-1) = sin 1/n wt (0) = sin nwt + sin 1/n wt : ');
if length(iWhichFFT) == 0
  iWhichFFT = 10;
end

dT = 1/12;        T = (1 : iMax)*dT;
dF = 2*pi/max(T); F = (1 : iMax)*dF;
dF = 1/max(T);    F = (1 : iMax); F = (F-1)*dF;
  
xx = 10*sin(2*pi*T) + 20*cos(2*pi*2*T) + 5*cos(2*pi/2*T) + rand(size(T)); plot(T,xx)
yy = fft(xx); plot(F,abs(yy))

while length(intersect(iWhichFFT,[10 1 -1 0])) == 1
  if iWhichFFT == +10
    %% no detrending
    fftspectrum_clr    = wfftspectrum_clr;
    fftspectrum_cld    = wfftspectrum_cld;
    fftspectrum_stemp  = wfftspectrum_stemp;
    fftspectrum_RHSurf = wfftspectrum_RHSurf;
  elseif iWhichFFT == 1
    fftspectrum_clr    = xfftspectrum_clr;
    fftspectrum_cld    = xfftspectrum_cld;
    fftspectrum_stemp  = xfftspectrum_stemp;
    fftspectrum_RHSurf = xfftspectrum_RHSurf;
  elseif iWhichFFT == -1
    fftspectrum_clr    = yfftspectrum_clr;
    fftspectrum_cld    = yfftspectrum_cld;
    fftspectrum_stemp  = yfftspectrum_stemp;
    fftspectrum_RHSurf = yfftspectrum_RHSurf;
  elseif iWhichFFT == 0
    fftspectrum_clr    = zfftspectrum_clr;
    fftspectrum_cld    = zfftspectrum_cld;
    fftspectrum_stemp  = zfftspectrum_stemp;
    fftspectrum_RHSurf = zfftspectrum_RHSurf;
  end
  
  mag_fftspectrum_clr = abs(fftspectrum_clr);
  phase_fftspectrum_clr = atan(imag(fftspectrum_clr)./real(fftspectrum_clr))*180/pi;
  mag_fftspectrum_cld = abs(fftspectrum_cld);
  phase_fftspectrum_cld = atan(imag(fftspectrum_cld)./real(fftspectrum_cld))*180/pi;
  
  mag_fftspectrum_stemp = abs(fftspectrum_stemp);
  phase_fftspectrum_stemp = atan(imag(fftspectrum_stemp)./real(fftspectrum_stemp))*180/pi;
  mag_fftspectrum_RHSurf = abs(fftspectrum_RHSurf);
  phase_fftspectrum_RHSurf = atan(imag(fftspectrum_RHSurf)./real(fftspectrum_RHSurf))*180/pi;
  
  figure(1); pcolor(log10(mag_fftspectrum_clr(:,(2:iMax)))); colorbar; shading flat;colormap jet; title('magnitude FFT');
  figure(2); pcolor(phase_fftspectrum_clr(:,1:iMax)); colorbar; shading flat;colormap jet; title('Phase FFT');

  figure(3); plot(F(2:iMax),log10(mag_fftspectrum_clr(:,(2:iMax))'),'c',F(2:iMax),log10(nanmean(mag_fftspectrum_clr(:,(2:iMax)),1)'),'b'); colorbar; shading flat;colormap jet; axis([0 iMax/2 0 2])
  figure(4); plot(F(2:iMax),(phase_fftspectrum_clr(:,(2:iMax))'),'c',F(2:iMax),(nanmean(phase_fftspectrum_clr(:,(2:iMax)),1)'),'b'); colorbar; shading flat;colormap jet; axis([1 iMax/2 -90 +90])
  
  figure(3); plot(F(2:iMax),log10(nanmean(mag_fftspectrum_clr(:,2:iMax),1)'),'b',F(2:iMax),log10(nanmean(mag_fftspectrum_cld(:,2:iMax),1)'),'r',...
                  F(2:iMax),log10(nanmean(mag_fftspectrum_stemp(:,2:iMax),1)'),'k',F(2:iMax),log10(nanmean(mag_fftspectrum_RHSurf(:,2:iMax),1)'),'g'); axis([1 iMax/2 0 2]); title('Magnitude FFT')
  
  figure(4); plot(F(2:iMax),nanmean(phase_fftspectrum_clr(:,2:iMax),1)','b',F(2:iMax),nanmean(phase_fftspectrum_cld(:,2:iMax),1)','r',...
                  F(2:iMax),nanmean(phase_fftspectrum_stemp(:,2:iMax),1)','k',F(2:iMax),nanmean(phase_fftspectrum_RHSurf(:,2:iMax),1)','g'); axis([1 iMax/2 -50 +50]); title('Phase FFT');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  it0 = 2:iMax/2;
  it = find( F(it0) >= 0.5 & F(it0) <= 5);
  figure(3); clf; loglog(F(it0),nanmean(mag_fftspectrum_clr(:,it0),1),F(it0(it)),nanmean(mag_fftspectrum_clr(:,it0(it)),1),'bo'); 
  P1 = polyfit(log(F(it0(it))),log(nanmean(mag_fftspectrum_clr(:,it0(it)),1)),1); Q1 = polyval(P1,log(F(it0(it)))); Q1 = exp(Q1); hold on; loglog(F(it0(it)),Q1,'b'); hold off

  it0 = 2:iMax/2;
  it = find( F(it0) >= 0.5 & F(it0) <= 5);
  P1 = polyfit(log(F(it0(it))),log(nanmean(mag_fftspectrum_clr(:,it0(it)),1)),1); Q1 = polyval(P1,log(F(it0(it)))); Q1 = exp(Q1); 
  P2 = polyfit(log(F(it0(it))),log(nanmean(mag_fftspectrum_cld(:,it0(it)),1)),1); Q2 = polyval(P2,log(F(it0(it)))); Q2 = exp(Q2); 
  P3 = polyfit(log(F(it0(it))),log(nanmean(mag_fftspectrum_stemp(:,it0(it)),1)),1); Q3 = polyval(P3,log(F(it0(it)))); Q3 = exp(Q3); 
  P4 = polyfit(log(F(it0(it))),log(nanmean(mag_fftspectrum_RHSurf(:,it0(it)),1)),1); Q4 = polyval(P4,log(F(it0(it)))); Q4 = exp(Q4); 

  figure(3); clf; loglog(F(2:iMax/2),nanmean(mag_fftspectrum_clr(:,2:iMax/2),1),'b',F(2:iMax/2),nanmean(mag_fftspectrum_cld(:,2:iMax/2),1),'r',...
                         F(2:iMax/2),nanmean(mag_fftspectrum_stemp(:,2:iMax/2),1),'k',F(2:iMax/2),nanmean(mag_fftspectrum_RHSurf(:,2:iMax/2),1),'g','linewidth',2); axis([0 10  0 100])
    hold on; loglog(F(it0(it)),Q1,'b',F(it0(it)),Q2,'r',F(it0(it)),Q3,'k',F(it0(it)),Q4,'g'); hold off; title('Magnitude FFT')
    hl = legend('CLR OLR','CLD OLR','STEMP','RH Surf','location','best','fontsize',10);
    str1 = ['CLR OLR ' num2str(P1(1))];
    str2 = ['CLD OLR ' num2str(P2(1))];
    str3 = ['STEMP '   num2str(P3(1))];
    str4 = ['RHSurf '  num2str(P4(1))];
    hl = legend(str1,str2,str3,str4,'location','best','fontsize',10);
    xlabel('frequency 1/yr'); ylabel('abs(FFT)')
 
  %%%%%%%%%%%%%%%%%%%%%%%%%

  iWhichFFT = input('Enter which fft (+10,default) = no detrend (1) = use sin nwt (-1) = sin 1/n wt (0) = sin nwt + sin 1/n wt : ');
  if length(iWhichFFT) == 0
    iWhichFFT = 10;
  end

end
