addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME

fin{1} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2012_clear_desc_ocean.mat';
fin{2} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2013_clear_desc_ocean.mat';
fin{3} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2014_clear_desc_ocean.mat';
fin{4} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2015_clear_desc_ocean.mat';
fin{5} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2016_clear_desc_ocean.mat';
fin{6} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2017_clear_desc_ocean.mat';
fin{7} = '/asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_2018_clear_desc_ocean.mat';
fin{8} =  '/home/sbuczko1/Work/cris_stats/rtp_cris_lowres_rad_2019_clear_desc_ocean.mat';
fin{9} =  '/home/sbuczko1/Work/cris_stats/rtp_cris_lowres_rad_2020_clear_desc_ocean.mat';
fin{10} = '/home/sbuczko1/Work/cris_stats/rtp_cris_lowres_rad_2021_clear_desc_ocean.mat';

%{
[sergio@taki-usr2 CRIS_new_clear_scan_January2020]$ ls -lt  /home/sbuczko1/Work/cris_stat*
total 9096060
-rw-rw---- 1 sbuczko1 pi_strow 2291186761 May 17 08:46 rtp_cris_lowres_rad_2020_clear_asc_ocean.mat
-rw-rw---- 1 sbuczko1 pi_strow 1788626872 May 17 08:37 rtp_cris_lowres_rad_2019_clear_asc_ocean.mat
-rw-rw---- 1 sbuczko1 pi_strow 2271897529 May 17 07:52 rtp_cris_lowres_rad_2020_clear_desc_ocean.mat
-rw-rw---- 1 sbuczko1 pi_strow 1779659266 May 17 07:46 rtp_cris_lowres_rad_2019_clear_desc_ocean.mat
-rw-rw---- 1 sbuczko1 pi_strow  592882053 May 17 06:27 rtp_cris_lowres_rad_2021_clear_asc_ocean.mat
-rw-rw---- 1 sbuczko1 pi_strow  590099224 May 17 06:07 rtp_cris_lowres_rad_2021_clear_desc_ocean.mat
%}

f = instr_chans('cris1317');
i1231 = find(f >= 1231,1);
i900 = find(f >= 900,1);
i1419 = find(f >= 1419,1);

%robs = [];
%rclr = [];
%rtime = [];
for ii = 1:10
  fprintf(1,'%s \n',fin{ii});
  a = load(fin{ii});
  wah = size(a.robs);
  if wah(4) == 1317
    ndays(ii) = wah(1);
    robs = [robs; squeeze(a.robs(:,:,5,[i900 i1231 i1419]))];
    if isfield(a,'rclr')
      rclr = [rclr; squeeze(a.rclr(:,:,5,[i900 i1231 i1419]))];
    else
      rclr = [rclr; 0.25*squeeze(a.robs(:,:,5,[i900 i1231 i1419]))];
    end
    rtime = [rtime; a.rtime_mean(:,:,5)];
  elseif wah(3) == 1317
    ndays(ii) = wah(1);
    robs = [robs; squeeze(a.robs(:,:,[i900 i1231 i1419],5))];
    if isfield(a,'rcal')
      rclr = [rclr; squeeze(a.rcal(:,:,[i900 i1231 i1419],5))];
    else
      rclr = [rclr; 0.25*squeeze(a.robs(:,:,[i900 i1231 i1419],5))];
    end
    rtime = [rtime; a.rtime_mean(:,:,5)];
  else
    error('oops')
  end
  %disp('ret to continue'); pause
  %whos robs rclr rtime
end

 [(2012:2021); ndays]'

whos robs rclr rtime
[yy,mm,dd,hh] = tai2utcSergio(mean(rtime,2));

yyx = yy + (mm-1)/12 + (dd-1)/30/12;
figure(1); plot(yyx,smooth(rad2bt(900,squeeze(robs(:,20,1))),90),yyx,smooth(rad2bt(900,squeeze(rclr(:,20,1))),90)); ylim([290 300])

