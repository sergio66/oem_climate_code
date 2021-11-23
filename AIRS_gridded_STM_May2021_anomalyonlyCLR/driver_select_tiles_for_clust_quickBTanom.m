addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /home/sergio/MATLABCODE/COLORMAP

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
ctmp = load('/home/sergio/MATLABCODE/PLOTTER/coast.mat');

iVers = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iVers == 1
  %% india
  oo1 = find(p.rlat >= +5 & p.rlat <= +35 & p.rlon >= +70 & p.rlon <= +90);
  
  %% Nino SST Region 3.4
  %% https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni
  oo2 = find(p.rlat >= -5 & p.rlat <= +5 & p.rlon >= -170 & p.rlon <= -120);
  
  oo = [oo1 oo2];
  
  whos oo
  plot(p.rlon,p.rlat,'.',p.rlon(oo),p.rlat(oo),'rx')
  hold on; plot(ctmp.long, ctmp.lat, 'k','linewidth',2); hold off
  
  fid = fopen('do_these_tiles.txt1','w');
  fprintf(fid,'%5i \n',oo);
  fclose(fid);

elseif iVers == 2
  %% artic oscillation
  oo1 = find(p.rlat >= -90 & p.rlat <= +90 & p.rlon >= -60 & p.rlon <= -20);
  
  %% southern ocean
  oo2 = find(p.rlat >= -75 & p.rlat <= -55 & p.rlon >= -180 & p.rlon <= +180);
  
  oo = [oo1 oo2];
  
  whos oo
  plot(p.rlon,p.rlat,'.',p.rlon(oo),p.rlat(oo),'rx')
  hold on; plot(ctmp.long, ctmp.lat, 'k','linewidth',2); hold off
  
  fid = fopen('do_these_tiles2.txt','w');
  fprintf(fid,'%5i \n',oo);
  fclose(fid);

elseif iVers == 3
  %% wv taperecorder around equator
  oo2 = find(p.rlat >= -8 & p.rlat <= +8 & p.rlon >= -180 & p.rlon <= +180);
  
  oo = [oo2];
  
  whos oo
  plot(p.rlon,p.rlat,'.',p.rlon(oo),p.rlat(oo),'rx')
  hold on; plot(ctmp.long, ctmp.lat, 'k','linewidth',2); hold off
  
  fid = fopen('do_these_tiles3.txt','w');
  fprintf(fid,'%5i \n',oo);
  fclose(fid);
end
