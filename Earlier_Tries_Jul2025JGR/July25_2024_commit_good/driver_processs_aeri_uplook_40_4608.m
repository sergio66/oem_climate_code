[h,ha,p40,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm_uplook.rtp');

iMonthSoFar = 120;
if ~exist('pnew_op')
  junk = load(['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(iMonthSoFar) '.mat'],'pnew_op');
  pnew_op = junk.pnew_op;
  clear junk
end

umbc_ecrad40 = load('ilrcalcs_CKD32.mat');
ecrad_4608 = load('ilr_4608_ecRad.mat');

kc40_0 = umbc_ecrad40.ilr0;
ind = (1:40) + 0*40; kc40_WV    = umbc_ecrad40.ilrX(ind) - kc40_0;
ind = (1:40) + 1*40; kc40_CO2   = umbc_ecrad40.ilrX(ind) - kc40_0;
ind = (1:40) + 2*40; kc40_T     = umbc_ecrad40.ilrX(ind) - kc40_0;
ind = (1:40) + 3*40; kc40_total = umbc_ecrad40.ilrX(ind) - kc40_0;

ecRad40_0 = umbc_ecrad40.ecRad0.clr;
ind = (1:40) + 0*40; ecRad40_WV    = umbc_ecrad40.ecRad.clr(ind) - ecRad40_0;
ind = (1:40) + 1*40; ecRad40_CO2   = umbc_ecrad40.ecRad.clr(ind) - ecRad40_0;
ind = (1:40) + 2*40; ecRad40_T     = umbc_ecrad40.ecRad.clr(ind) - ecRad40_0;
ind = (1:40) + 3*40; ecRad40_total = umbc_ecrad40.ecRad.clr(ind) - ecRad40_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); plot(p40.rlat,kc40_WV,'b',p40.rlat,kc40_CO2,'g',p40.rlat,kc40_T,'r',p40.rlat,kc40_total,'k','linewidth',2)
  plotaxis2; hl = legend('WV','CO2','T','total','location','best','fontsize',10); 
  xlabel('Latitude'); ylabel('W/m2/yr'); title('kCARTA 40 profiles')
  ylim([0 0.2])

figure(2); plot(p40.rlat,ecRad40_WV,'b',p40.rlat,ecRad40_CO2,'g',p40.rlat,ecRad40_T,'r',p40.rlat,ecRad40_total,'k','linewidth',2)
  plotaxis2; hl = legend('WV','CO2','T','total','location','best','fontsize',10); 
  xlabel('Latitude'); ylabel('W/m2/yr'); title('ecRad 40 profiles')
  ylim([0 0.2])

figure(3); plot(pnew_op.rlat,ecrad_4608.ecRadWV,'b',pnew_op.rlat,ecrad_4608.ecRadCO2,'g',pnew_op.rlat,ecrad_4608.ecRadT,'r',pnew_op.rlat,ecrad_4608.ecRadtotal,'k','linewidth',2)
  plotaxis2; hl = legend('WV','CO2','T','total','location','best','fontsize',10); 
  xlabel('Latitude'); ylabel('W/m2/yr'); title('ecRad 40 profiles')
  ylim([0 0.2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cf = [0.10 0.20];  %% coeffs for quadratic fit
cf = [0.05 0.10];  %% coeffs for quadratic fit
cf = [0.06 0.17];  %% coeffs for quadratic fit
figure(4); plot(p40.rlat,kc40_WV./kc40_total,'b',p40.rlat,kc40_CO2./kc40_total,'g',p40.rlat,kc40_T./kc40_total,'r',p40.rlat,kc40_total./kc40_total,'k',p40.rlat,cf(1) + cf(2)*sin(pi/180*p40.rlat).^2,'g--','linewidth',2)
  plotaxis2; hl = legend('WV','CO2','T','total','location','best','fontsize',10); 
  xlabel('Latitude'); ylabel('W/m2/yr'); title('kCARTA 40 profiles')
  ylim([0 1])

figure(5); plot(p40.rlat,ecRad40_WV./ecRad40_total,'b',p40.rlat,ecRad40_CO2./ecRad40_total,'g',p40.rlat,ecRad40_T./ecRad40_total,'r',p40.rlat,ecRad40_total./ecRad40_total,'k',...
                p40.rlat,cf(1) + cf(2)*sin(pi/180*p40.rlat).^2,'g--','linewidth',2)
  plotaxis2; hl = legend('WV','CO2','T','total','location','best','fontsize',10); 
  xlabel('Latitude'); ylabel('W/m2/yr'); title('ecRad 40 profiles')
  ylim([0 1])

ecRadtotal = ecrad_4608.ecRadtotal;
slat4608 = sin(pnew_op.rlat*pi/180);
co24608  = ecrad_4608.ecRadCO2./ecRadtotal;
n = input('Enter order of polynomial fit : ');
P = polyfit(slat4608,co24608,n)
figure(6); plot(pnew_op.rlat,ecrad_4608.ecRadWV./ecRadtotal,'b',pnew_op.rlat,ecrad_4608.ecRadCO2./ecRadtotal,'c',pnew_op.rlat,ecrad_4608.ecRadT./ecRadtotal,'r',pnew_op.rlat,ecrad_4608.ecRadtotal./ecRadtotal,'k',...
                pnew_op.rlat,cf(1) + cf(2)*sin(pi/180*pnew_op.rlat).^2,'g--',pnew_op.rlat,polyval(P,slat4608),'gx-','linewidth',2)
  plotaxis2; hl = legend('WV','CO2','T','total','location','best','fontsize',10); 
  xlabel('Latitude'); ylabel('W/m2/yr'); title('ecRad 40 profiles')
  ylim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% >> driver_processs_aeri_uplook_40_4608
% Enter order of polynomial fit : >>3
% P = -0.0149    0.1771   -0.0486    0.0562
% >> driver_processs_aeri_uplook_40_4608
% Enter order of polynomial fit : >>2
% P =    0.1771   -0.0595    0.0562
