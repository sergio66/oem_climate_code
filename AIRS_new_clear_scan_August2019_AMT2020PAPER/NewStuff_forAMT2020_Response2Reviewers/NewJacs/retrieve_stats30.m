addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE

clear all

sarta   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

iType = 1;  %% works, very safe
iType = 2;  %% works???

if iType == 1

  stats30 = load('/home/strow/Work/Airs/Stability/Data_old/Desc/statlat30.mat');
  %{
  count       5065x2378            96356560  double
  gas1        5065x60               1215600  single
  gas3        5065x60               1215600  single
  iudef4      5065x1                  40520  double
  lat         5065x1                  20260  single
  lon         5065x1                  20260  single
  nlevs       5065x1                  40520  double
  plevs       5065x60               1215600  single
  ptemp       5065x60               1215600  single
  rcal        5065x2378            48178280  single
  robs        5065x2378            48178280  single
  rtime       5065x1                  40520  double
  satzen      5065x1                  20260  single
  solzen      5065x1                  20260  single
  spres       5065x1                  20260  single
  stemp       5065x1                  20260  single
  %}

  good = find(stats30.nlevs > 0);

  [h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');

  [h,p] = replicate_rtp_headprof(h,pbefore,30,length(good));  

  p.gas_1 = stats30.gas1(good,:)';
  p.gas_3 = stats30.gas3(good,:)';
  p.ptemp = stats30.ptemp(good,:)';
  p.plevs = stats30.plevs(good,:)';
  p.nlevs = int32(stats30.nlevs(good)');
  p.spres = stats30.spres(good)';
  p.stemp = stats30.stemp(good)';
  p.satzen = stats30.satzen(good)';
  p.solzen = stats30.solzen(good)';
  p.robs1 = stats30.robs(good,:)';
  p.rcalc = stats30.rcal(good,:)';
  p.rtime = stats30.rtime(good)';

  p = rmfield(p,'gas_2');
  p = rmfield(p,'gas_4');
  p = rmfield(p,'gas_5');
  p = rmfield(p,'gas_6');
  p = rmfield(p,'gas_9');
  p = rmfield(p,'gas_12');
  p = rmfield(p,'palts');

  h.ptype = 0;
  h.ngas = 2;
  h.gunit = [21 21]';
  h.glist = [ 1  3]';
  h.pfields = 7;
  rtpwrite('stats_latbin30.ip.rtp',h,ha,p,pa);

elseif iType == 2
  stats30 = load('/home/strow/Work/Airs/Stability/Data/Desc/statlat30.mat');
  %{
           count: [6152x2645 double]
       gas1_mean: [6152x101 double]
       gas3_mean: [6152x101 double]
     iudef4_mean: [6152x1 double]
        lat_mean: [6152x1 double]
        lon_mean: [6152x1 double]
      nlevs_mean: [6152x1 double]
      ptemp_mean: [6152x101 double]
            rclr: [6152x2645 double]
    rclrbias_std: [6152x2645 double]
            robs: [6152x2645 double]
      rtime_mean: [6152x1 double]
     satzen_mean: [6152x1 double]
     solzen_mean: [6152x1 double]
      spres_mean: [6152x1 double]
      stemp_mean: [6152x1 double]
        tcc_mean: [6152x1 double]
  %}

  good = find(stats30.nlevs_mean > 0);

%{
  [h,ha,pbefore,pa] = rtpread('/asl/rtp/rtp_airicrad_v6/clear/2002/era_airicrad_day244_clear.rtp');
  boo = find(pbefore.solzen > 90 & pbefore.landfrac == 0 & abs(pbefore.rlat) <= 20,1);
  [h,pbefore] = subset_rtp(h,pbefore,[],[],boo);
  rtpwrite('junk.ip.rtp',h,ha,pbefore,pa);
  klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp'];
  eval(klayerser);
  [h,ha,pbefore,pa] = rtpread('junk.op.rtp');
%}

  [h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');
  pbefore.solzen = 150 * ones(size(pbefore.stemp));
  boo = find(pbefore.solzen > 90 & pbefore.landfrac == 0 & abs(pbefore.rlat) <= 20,1);  
  [h,pbefore] = subset_rtp(h,pbefore,[],[],boo);

  [h,p] = replicate_rtp_headprof(h,pbefore,1,length(good));  

  p.gas_1 = stats30.gas1_mean(good,:)';
  p.gas_3 = stats30.gas3_mean(good,:)';
  p.ptemp = stats30.ptemp_mean(good,:)';
  p.nlevs = int32(stats30.nlevs_mean(good)');
  p.spres = stats30.spres_mean(good)';
  p.stemp = stats30.stemp_mean(good)';
  p.satzen = stats30.satzen_mean(good)';
  p.solzen = stats30.solzen_mean(good)';
  p.robs1 = stats30.robs(good,:)';
  p.rcalc = stats30.rclr(good,:)';
  p.rtime = stats30.rtime_mean(good)';

  boo = find(isnan(p.ptemp) | isnan(p.gas_1) | isnan(p.gas_3));
  p.ptemp(boo) = 300;
  p.gas_1(boo) = 0.0;
  p.gas_3(boo) = 0.0;

  [hjunk,ha,pjunk,pa] = rtpread('/asl/rtp/rtp_airicrad_v6/clear/2002/era_airicrad_day244_clear.rtp');  
  h.nchan = hjunk.nchan;
  h.ichan = hjunk.ichan;
  h.vchan = hjunk.vchan;

  rtpwrite('stats_latbin30.op.rtp',h,ha,p,pa);
end

klayerser = ['!' klayers ' fin=stats_latbin30.ip.rtp fout=stats_latbin30.op.rtp'];
sartaer =   ['!' sarta ' fin=stats_latbin30.op.rtp fout=stats_latbin30.rp.rtp'];
if iType == 1
  eval(klayerser);
end
eval(sartaer);

[hnew,ha,pnew,pa] = rtpread('stats_latbin30.rp.rtp');
tobs = rad2bt(h.vchan,p.robs1);
tcal0 = rad2bt(h.vchan,p.rcalc);
tcalF = rad2bt(h.vchan,pnew.rcalc);

plot(hnew.vchan,mean(tobs'-tcal0'),'b',hnew.vchan,mean(tobs'-tcalF'),'r',...
    hnew.vchan,std(tobs'-tcal0'),'c--',hnew.vchan,std(tobs'-tcalF'),'m--')

iCut = 20;
iNum = ceil(length(pnew.rtime)/iCut);
for ii = 1 : 20
  ind = (1:iNum) + (ii-1)*iNum;
  yes = find(ind <= length(pnew.rtime));
  ind = ind(yes);
  fprintf(1,'%3i : %4i to %4i \n',ii,ind(1),ind(end));
  [hsavei,psavei] = subset_rtp(hnew,pnew,[],[],ind);
  fout = ['CUTUP_latbin30/cutup_latbin30_' num2str(ii) '.rtp'];
  rtpwrite(fout,hsavei,ha,psavei,pa);
end
