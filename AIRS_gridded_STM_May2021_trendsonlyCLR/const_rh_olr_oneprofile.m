%function olr = const_rh_olr;

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio//MATLABCODE/
addpath /home/sergio//MATLABCODE/CONVERT_GAS_UNITS/
addpath /home/sergio//MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/

disp('do this before starting Matlab : module load netCDF-Fortran/4.4.4-intel-2018b')

fip = 'const_rh_olr_oneprofile.ip.rtp';
fop = 'const_rh_olr_oneprofile.op.rtp';
frp = 'const_rh_olr_oneprofile.rp.rtp';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
rmer    = ['!/bin/rm ' fip ' ' fop ' ' frp]; 

fin = '/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_g80_op.so2.rtp_new';
fin = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm.rtp';
[h,ha,p,pa] = rtpread(fin);

frac = 2.0; % increase water by too much   to mantain const RH
frac = 0.5; % increase water by too little to mantain const RH

h.vchan = instr_chans;
[h,p] = subset_rtp(h,p,[],[],20);
p.salti = 0;
p.landfrac = 0;
%p.ptemp = p.ptemp(1:98);
%p.gas_1 = p.gas_1(1:98);
%p.gas_2 = p.gas_2(1:98);
%p.gas_3 = p.gas_3(1:98);
%p.gas_4 = p.gas_4(1:98);
%p.gas_5 = p.gas_5(1:98);
%p.gas_6 = p.gas_6(1:98);
%p.gas_9 = p.gas_9(1:98);
%p.gas_12 = p.gas_12(1:98);

dT = 1.0;
dT = 0.025;

figure(1)
ones101 = linspace(-1,+1,101)';
ones101 = dT*ones101;
semilogy(ones101,p.plevs); plotaxis2; set(gca,'ydir','reverse'); ylim([0.01 p.spres])

i100 = find(p.plevs >= 100,1);
ones101 = [dT*(linspace(-1,+1,i100)-0.0) dT*ones(1,length([i100+1:101]))]';
ones101 = [dT*(linspace(-1,+1,i100)-0.5) dT*ones(1,length([i100+1:101]))]';
ones101 = [dT*(linspace(-1,+1,i100)-0.5) dT*linspace(-0.2,+1,101-i100)]';

slope = (dT*0.5-dT)/(101-i100); intercept = -dT*0.5 - slope*1; rest = (1:101-i100)*slope + intercept;
ones101 = [dT*(linspace(-1,+1,i100)-0.5) -rest]';
semilogy(ones101,p.plevs); plotaxis2; set(gca,'ydir','reverse'); ylim([0.01 p.spres])

%% initial state, TROPICAL
px = p;
rh0 = layeramt2RH(h,px,1);
olr0 = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olr0.rcalc = pjunk.rcalc;

%% perturb surf temp and air temp
px = p;
px.ptemp = px.ptemp + ones101;
px.stemp = px.stemp + dT;
px.gas_2 = px.gas_2 * (1+2.2/385);
rhTA = layeramt2RH(h,px,1);
olrTA = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olrTA.rcalc = pjunk.rcalc;

%% perturb surf temp only
px = p;
px.stemp = px.stemp + dT;
px.gas_2 = px.gas_2 * (1+2.2/385);
rhTS = layeramt2RH(h,px,1);
olrTS = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olrTS.rcalc = pjunk.rcalc;

%% purely water vapor only
px = p;
px.gas_2 = px.gas_2 * (1+2.2/385);
wvMult = 1-(rhTA-rh0)./rh0; wvMult(101) = 1;
px.gas_1 = px.gas_1 .* wvMult;
rhWV0 = layeramt2RH(h,px,1);
olrWV0 = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olrWV0.rcalc = pjunk.rcalc;

%% purely half water vapor only
px = p;
px.gas_2 = px.gas_2 * (1+2.2/385);
wvMult = 1-(rhTA-rh0)./rh0*frac; wvMult(101) = 1;
px.gas_1 = px.gas_1 .* wvMult;
rhWV0x = layeramt2RH(h,px,1);
olrWV0x = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olrWV0x.rcalc = pjunk.rcalc;

%% perturb surf temp and air temp and WV to get d(RH) ~ 0
px = p;
px.gas_2 = px.gas_2 * (1+2.2/385);
px.ptemp = px.ptemp + ones101;
px.stemp = px.stemp + dT;
wvMult = 1-(rhTA-rh0)./rh0; wvMult(101) = 1;
px.gas_1 = px.gas_1 .* wvMult;
rhWV = layeramt2RH(h,px,1);
olrWV = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olrWV.rcalc = pjunk.rcalc;

%% perturb surf temp and air temp and WV by half
px = p;
px.gas_2 = px.gas_2 * (1+2.2/385);
px.ptemp = px.ptemp + ones101;
px.stemp = px.stemp + dT;
wvMult = 1-(rhTA-rh0)./rh0 * frac; wvMult(101) = 1;
px.gas_1 = px.gas_1 .* wvMult;
rhWVx = layeramt2RH(h,px,1);
olrWVx = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
rtpwrite(fop,h,[],px,[]);
eval(sartaer);
[hjunk,~,pjunk,~] = rtpread(frp);
olrWVx.rcalc = pjunk.rcalc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); 
playsN = p.plevs(1:100)-p.plevs(2:101);
playsD = log(p.plevs(1:100)./p.plevs(2:101));
p.plays = playsN./playsD;
semilogy(rhTS-rh0,p.plays,rhTA-rh0,p.plays,rhWV-rh0,p.plays,'x-',rhWVx-rh0,p.plays,'linewidth',2); set(gca,'ydir','reverse'); ylim([10 1000]);
  hl = legend('TS only','TS+T(z)','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  ylabel('P(mb)'); xlabel('\delta RH = RHx - RH0')

figure(2);
bands = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
bands17 = bands(1:end-1)+bands(2:end); bands17 = bands17/2;
plot(bands17,olr0.bands.clr,bands17,olrTS.bands.clr,bands17,olrTA.bands.clr,bands17,olrWV.bands.clr,'linewidth',2)
plot(bands17,olrTS.bands.clr-olr0.bands.clr,bands17,olrTA.bands.clr-olr0.bands.clr,bands17,olrWV.bands.clr-olr0.bands.clr,'x-',bands17,olrWVx.bands.clr-olr0.bands.clr,'linewidth',2)
  hl = legend('TS only','TS+T(z)','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  xlabel('RRTM Center Wavenumber cm^{-1}'); ylabel('\delta Flux = Fx - F0 (W/m2)');

eval(rmer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); 
playsN = p.plevs(1:100)-p.plevs(2:101);
playsD = log(p.plevs(1:100)./p.plevs(2:101));
p.plays = playsN./playsD;
semilogy(rhTS-rh0,p.plays,rhTA-rh0,p.plays,rhWV0-rh0,p.plays,rhWV-rh0,p.plays,'x-',rhWVx-rh0,p.plays,'linewidth',2); set(gca,'ydir','reverse'); ylim([10 1000]);
  hl = legend('TS only','TS+T(z)','increase WV only','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  ylabel('P(mb)'); xlabel('\delta RH = RHx - RH0')

figure(2);
bands = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
bands17 = bands(1:end-1)+bands(2:end); bands17 = bands17/2;
plot(bands17,olr0.bands.clr,bands17,olrTS.bands.clr,bands17,olrTA.bands.clr,bands17,olrWV.bands.clr,'linewidth',2)
plot(bands17,olrTS.bands.clr-olr0.bands.clr,bands17,olrTA.bands.clr-olr0.bands.clr,bands17,olrWV0.bands.clr-olr0.bands.clr,bands17,olrWV.bands.clr-olr0.bands.clr,'x-',bands17,olrWVx.bands.clr-olr0.bands.clr,'linewidth',2)
  hl = legend('TS only','TS+T(z)','increase WV only','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  xlabel('RRTM Center Wavenumber cm^{-1}'); ylabel('\delta Flux = Fx - F0 (W/m2)');

BT(:,1) = rad2bt(hjunk.vchan,olr0.rcalc);
BT(:,2) = rad2bt(hjunk.vchan,olrTS.rcalc);
BT(:,3) = rad2bt(hjunk.vchan,olrTA.rcalc);
BT(:,4) = rad2bt(hjunk.vchan,olrWV0.rcalc);
BT(:,5) = rad2bt(hjunk.vchan,olrWV0x.rcalc); %% not really used
BT(:,6) = rad2bt(hjunk.vchan,olrWV.rcalc);
BT(:,7) = rad2bt(hjunk.vchan,olrWVx.rcalc);
figure(3);
plot(hjunk.vchan,BT(:,[2 3 4 6 7]) - BT(:,1)*ones(1,5),'linewidth',2);
hold on; plot(hjunk.vchan,BT(:,6)-BT(:,1),'x','color',[0.4940, 0.1840, 0.5560]); hold off
  plotaxis2; xlim([640 1640]); hl = legend('TS only','TS+T(z)','increase WV only','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  xlabel('AIRS Wavenumber cm^{-1}'); ylabel('\delta BT = BTx - BT0 (K)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf
ta = tiledlayout(1,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
semilogy(rhTS-rh0,p.plays,rhTA-rh0,p.plays,rhWV0-rh0,p.plays,rhWV-rh0,p.plays,'x-',rhWVx-rh0,p.plays,'linewidth',2); set(gca,'ydir','reverse'); ylim([10 1000]);
  %hl = legend('TS only','TS+T(z)','increase WV only','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  ylabel('P(mb)'); xlabel('\delta RH = RHx - RH0')

tafov(2) = nexttile;
plot(bands17,olrTS.bands.clr-olr0.bands.clr,bands17,olrTA.bands.clr-olr0.bands.clr,bands17,olrWV0.bands.clr-olr0.bands.clr,bands17,olrWV.bands.clr-olr0.bands.clr,'x-',bands17,olrWVx.bands.clr-olr0.bands.clr,'linewidth',2)
  %hl = legend('TS only','TS+T(z)','increase WV only','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
  xlabel('RRTM Center Wavenumber cm^{-1}'); ylabel('\delta Flux = Fx - F0 (W/m2)');

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

%hl = legend('TS only','TS+T(z)','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',10);
%lg = legend('TS only','TS+T(z)','TS+T(z),increase WV \newline \delta RH \rightarrow 0','TS+T(z), incr WV by half','location','best','fontsize',8);
lg = legend('TS only','TS+T(z)','increase WV only','TS+T(z),increase WV','TS+T(z), increase WV by 1/2','location','best','fontsize',8);
lg.Layout.Tile = 'North';

figure(5);
plot(hjunk.vchan,BT(:,[6 7]) - BT(:,1)*ones(1,2),'linewidth',2);
  plotaxis2; xlim([640 1640]); hl = legend('RH constant','RH $\downarrow$','location','best','fontsize',10,'interpreter','latex');
  xlabel('Wavenumber cm^{-1}'); ylabel('\delta BT = BTx - BT0 (K)');

% https://www.mathworks.com/matlabcentral/answers/160487-how-can-i-draw-a-line-with-arrow-head-between-2-data-points-in-a-plot
ylim([-0.1 +0.075])
%% get info specific to the axes you plan to plot into
set(gcf,'Units','normalized');
set(gca,'Units','normalized');
ax = axis;
ap = get(gca,'Position');

%% annotation from 640 to 800,0.05 to 0.05 1,2 to 3,4
xo = [640 800];
yo = [0.05 0.05];
xArrow = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
yArrow = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
%ah=annotation('textarrow',xArrow,yArrow,'String','T/CO2','FontSize',13,'Linewidth',2); ah.Color = 'red'; ah.Head1Style = 
ah=annotation('doublearrow',xArrow,yArrow); ah.Color = 'red'; text(670,0.06,'T/CO2','FontSize',13)

%% annotation from 1240 to 1640,0.05 to 0.05 1,2 to 3,4
xo = [1240 1640];
yo = [0.05 0.05];
xArrow = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
yArrow = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
ah=annotation('doublearrow',xArrow,yArrow); ah.Color = 'blue'; text(1400,0.06,'WV','FontSize',13)

%% annotation from 800 to 1000,0.05 to 0.05 1,2 to 3,4
xo = [800 1000];
yo = [0.05 0.05];
xArrow = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
yArrow = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
ah=annotation('doublearrow',xArrow,yArrow); ah.Color = 'green'; 
%% annotation from 800 to 1000,0.05 to 0.05 1,2 to 3,4
xo = [1100 1240];
yo = [0.05 0.05];
xArrow = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
yArrow = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
ah=annotation('doublearrow',xArrow,yArrow); ah.Color = 'green'; 
text(950,0.06,'WINDOW','FontSize',13)

%xArrow = [0  (800-640)/(1640-640)];  yArrow = [(0.05--0.1)/(0.075--0.1) (0.05--0.1)/(0.075--0.1)];annotation('textarrow',xArrow,yArrow,'String','T/CO2','FontSize',13,'Linewidth',2)
%text(800,0.05,'\xleftrightarrow','fontsize',20); text(800,0.05,'\xleftrightarrow','fontsize',20); 

addpath /asl/matlib/plotutils
%{
figure(3); aslprint('Figs_STMOct2021/delta_BT_with_WV.pdf')
figure(4); aslprint('Figs_STMOct2021/delta_OLR_with_WV.pdf')
figure(5); aslprint('Figs_STMOct2021/delta_BT_with_WVv2.pdf')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'lamba TS  = %8.6f \n',-(olrTS.clr-olr0.clr)/dT)
fprintf(1,'lamba TA  = %8.6f \n',-(olrTA.clr-olr0.clr)/dT)
%fprintf(1,'lamba WV  = %8.6f \n',-(olrWV.clr-olrTA.clr)/dT)
%fprintf(1,'lamba WVx = %8.6f \n',-(olrWVx.clr-olrTA.clr)/dT)
fprintf(1,'lamba WV  = %8.6f \n',-(olrWV0.clr-olrTA.clr)/dT)
fprintf(1,'lamba WVx = %8.6f \n',-(olrWV0x.clr-olrTA.clr)/dT)
