addpath /asl/matlib/h4tools

[h4608,~,p4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002.rtp');
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  ind = (1:72) + (ii-1)*72;
  [hjunk,pjunk] = subset_rtp_allcloudfields(h4608,p4608,[],[],ind);

  ix = find(pjunk.cfrac <= 0);
  pjunk.cfrac(ix) = NaN;
  pjunk.cprtop(ix) = NaN;
  pjunk.cprbot(ix) = NaN;
  pjunk.ctype(ix) = NaN;
  pjunk.cngwat(ix) = NaN;
  pjunk.cpsize(ix) = NaN;

  ix = find(pjunk.cfrac2 <= 0);
  pjunk.cfrac2(ix) = NaN;
  pjunk.cprtop2(ix) = NaN;
  pjunk.cprbot2(ix) = NaN;
  pjunk.ctype2(ix) = NaN;
  pjunk.cngwat2(ix) = NaN;
  pjunk.cpsize2(ix) = NaN;

  [pavgii] = find_average_rtp(pjunk);
  if ii == 1
    havg = hjunk;
    pavg = pavgii;
  else
    [havg,pavg] = cat_rtp(havg,pavg,hjunk,pavgii);
  end
end

figure(1); clf; plot(p4608.rlat,p4608.stemp,'b',pavg.rlat,pavg.stemp,'r'); title('stemp')
figure(2); clf; plot(p4608.rlat,p4608.cfrac,'b',pavg.rlat,pavg.cfrac,'c',p4608.rlat,p4608.cfrac2,'r',pavg.rlat,pavg.cfrac2,'m'); title('cfrac');
figure(3); clf; plot(p4608.rlat,p4608.cprtop,'b',pavg.rlat,pavg.cprtop,'c',p4608.rlat,p4608.cprtop2,'r',pavg.rlat,pavg.cprtop2,'m'); title('cprtop'); set(gca,'ydir','reverse');
figure(4); clf; plot(p4608.rlat,p4608.cngwat,'b',pavg.rlat,pavg.cngwat,'c',p4608.rlat,p4608.cngwat2,'r',pavg.rlat,pavg.cngwat2,'m'); title('cngwat');
figure(5); clf; plot(p4608.rlat,p4608.cpsize,'b',pavg.rlat,pavg.cpsize,'c',p4608.rlat,p4608.cpsize2,'r',pavg.rlat,pavg.cpsize2,'m'); title('cpsize');
figure(6); clf; plot(p4608.rlat,p4608.ctype,'b',pavg.rlat,pavg.ctype,'c',p4608.rlat,p4608.ctype2,'r',pavg.rlat,pavg.ctype2,'m'); title('ctype');

if ~exist('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_zonalavg64lats.rtp')
  rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_zonalavg64lats.rtp',havg,[],pavg,[]);
else
  disp('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_zonalavg64lats.rtp already exists')
end
