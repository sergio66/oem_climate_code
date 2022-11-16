clear all

addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools

%% now do quick sarta jacs
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;
fop = 'junk.op.rtp';
frp = 'junk.rp.rtp';

% BT1231cld quants = [0 0.03 0.05 0.1 0.2 0.5 0.8 0.9 0.95 0.97 1.0];
fop0 = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/2012/FixedNAN/all4608_era5_full12months_Qcumulative09.rtp';
[h0,ha,p0,pa] = rtpread(fop0);

iX = input('Enter [iLonBin iLatBin] : ');
iLonBin = 01; iLatBin = 32;
iLonBin = iX(1); iLatBin = iX(2);
miaow = (iLatBin-1)*72 + iLonBin; 
fprintf(1,'[iLonBin iLatBin miaow] = %2i %2i %4i \n',[iLonBin iLatBin miaow]);
fprintf(1,'[rLon rLat landfrac]    = %8.4f %8.4f %8.4f \n',[p0.rlon(miaow) p0.rlat(miaow) p0.landfrac(miaow)])
[h0,p0] = subset_rtp_allcloudfields(h0,p0,[],[],miaow);
hjunk = h0; 
[m_ts_jac,nlays,qrenorm,freq2645,colo3] = get_jac_fast([],miaow,iLonBin,iLatBin,2022);

pjunk = p0;
rtpwrite(fop,hjunk,[],pjunk,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp]; eval(sartaer);
[hjunk,~,pjunk0,~] = rtpread(frp);
bt0 = rad2bt(hjunk.vchan,pjunk0.rcalc);

pjunk = p0;
pjunk.gas_1 = pjunk.gas_1 * 1.1;
rtpwrite(fop,hjunk,[],pjunk,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp]; eval(sartaer);
[hjunk,~,pjunk1,~] = rtpread(frp);
bt1 = rad2bt(hjunk.vchan,pjunk1.rcalc);

pjunk = p0;
pjunk.gas_2 = pjunk.gas_2 * 1.1;
rtpwrite(fop,hjunk,[],pjunk,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp]; eval(sartaer);
[hjunk,~,pjunk2,~] = rtpread(frp);
bt2 = rad2bt(hjunk.vchan,pjunk2.rcalc);

pjunk = p0;
pjunk.ptemp = pjunk.ptemp + 0.1;
rtpwrite(fop,hjunk,[],pjunk,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp]; eval(sartaer);
[hjunk,~,pjunkT,~] = rtpread(frp);
btT = rad2bt(hjunk.vchan,pjunkT.rcalc);

pjunk = p0;
pjunk.stemp = pjunk.stemp + 0.1;
rtpwrite(fop,hjunk,[],pjunk,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp]; eval(sartaer);
[hjunk,~,pjunkT,~] = rtpread(frp);
btST = rad2bt(hjunk.vchan,pjunkT.rcalc);

%miaowx = 1;
%figure(1); clf; plot(hjunk.vchan,(btST(:,miaowx)-bt0(:,miaowx))/0.1);     xlim([640 1640]); title('ST jac')
%figure(1); clf; plot(hjunk.vchan,(btT(:,miaowx)-bt0(:,miaowx))/0.1);      xlim([640 1640]); title('T jac')
%figure(1); clf; plot(hjunk.vchan,(bt1(:,miaowx)-bt0(:,miaowx))/log(1.1)); xlim([640 1640]); title('WV jac')
%figure(1); clf; plot(hjunk.vchan,(bt2(:,miaowx)-bt0(:,miaowx))/log(1.1)); xlim([640 1640]); title('CO2 jac')

str = ['Lon/Lat/Prof = ' num2str(iLonBin,'%02d') ' ' num2str(iLatBin,'%02d') ' ' num2str(miaow,'%04d')];

pjunkT.robs1 = pjunkT.rcalc;
print_cloud_params(hjunk,pjunkT,1);

%% qrenrom(150) = 0.01 = T renorm;
ind = (6+1*97+(1:97)); figure(1); clf; plot(hjunk.vchan,0.01*(btT-bt0)/0.1,'b.-',hjunk.vchan,sum(m_ts_jac(:,ind),2),'r');       xlim([640 1640]); title([ str ' T jac'])
  hl = legend('SARTA','kCARTA','location','best');
ind = 6;               figure(2); clf; plot(hjunk.vchan,0.10*(btST-bt0)/0.1,'b.-',hjunk.vchan,sum(m_ts_jac(:,ind),2),'r');      xlim([640 1640]); title([ str ' ST jac'])
  hl = legend('SARTA','kCARTA','location','best');

%% qrenrom(90) = 0.01 = WV renorm;
ind = (6+0*97+(1:97)); figure(3); clf; plot(hjunk.vchan,0.01*(bt1-bt0)/log(1.1),'b.-',hjunk.vchan,sum(m_ts_jac(:,ind),2),'r'); xlim([640 1640]); title([ str ' WV jac'])
  hl = legend('SARTA','kCARTA','location','best');
ind = 1; figure(4); clf; plot(hjunk.vchan,2.2/400*(bt2-bt0)/log(1.1),'b.-',hjunk.vchan,m_ts_jac(:,ind),'r');                   xlim([640 1640]); title([ str ' CO2 jac'])
  hl = legend('SARTA','kCARTA','location','best');

ind = (6+2*97+(1:97)); figure(5); clf; plot(hjunk.vchan,0.01*colo3,'b.-',hjunk.vchan,sum(m_ts_jac(:,ind),2),'r');       xlim([640 1640]); title([ str ' O3 jac'])
  hl = legend('COLUMN','sum(O3jac(z))','location','best');
