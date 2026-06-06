%% see /home/sergio/git/matlabcode/find_average_rtp.m
%% first fix good old "all"   ptemp/gas_1,gas_3,cloud params as needed

pallx = pall;

[xmm,xnn,xoo] = size(pall.ptemp);
wonk = flipud(load('/home/sergio/git/matlabcode/airslevels.dat')) * ones(1,4608);
bonk = zeros(size(pall.ptemp));
for ii = 1 : xmm
  bonk(ii,:,:) = wonk;
end
pallx.plevs = bonk;
disp('made pallx.plevs')

pallx.nlevs = pallx.nlays + 1;

pallx.ptemp(pallx.ptemp <= 150) = NaN;
pallx.gas_1(pallx.gas_1 <=   0) = NaN;
pallx.gas_3(pallx.gas_3 <=   0) = NaN;

hallx.ngas = 2;
hallx.glist = [1 3];

disp('nan-ning plevs,ptemp,gas_N below p.nlevs')
[iLenZ,~] = size(pallx.plevs);
for ii = 1 : length(pallx.stemp)
  nlevs = pallx.nlevs(ii);
  pallx.plevs(nlevs+1:iLenZ,ii) = NaN;
  pallx.ptemp(nlevs+1:iLenZ,ii) = NaN;
  for gg = 1 : hallx.ngas
    str = ['pallx.gas_' num2str(hallx.glist(gg)) '(nlevs+1:iLenZ,ii) = NaN;'];
    eval(str)
  end
end

boo = find(pallx.ctype < 0 | pallx.cngwat < 0 | pallx.cpsize < 0 | pallx.cprtop < 0 | pallx.cprbot < 0);
  pallx.ctype(boo) = NaN;
  pallx.cngwat(boo) = NaN;
  pallx.cpsize(boo) = NaN;
  pallx.cprtop(boo) = NaN;
  pallx.cprbot(boo) = NaN;
boo = find(pallx.ctype < 0 | pallx.cngwat2 < 0 | pallx.cpsize2 < 0 | pallx.cprtop2 < 0 | pallx.cprbot2 < 0);
  pallx.ctype2(boo) = NaN;
  pallx.cngwat2(boo) = NaN;
  pallx.cpsize2(boo) = NaN;
  pallx.cprtop2(boo) = NaN;
  pallx.cprbot2(boo) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now average

xnames = fieldnames(pallx);
for ii = 1 : length(xnames)
  xxname = xnames{ii};
  str = ['junk = pallx.' xxname ';'];
  eval(str);
  [xmm,xnn,xoo] = size(junk);
  if xmm == 1 & xoo == 1
    xavg = nanmean(junk);
    xavg = junk;    
    str = ['avg_pall.' xxname ' = xavg;'];
    eval(str);
  elseif xmm >  1 & xoo == 1
    xavg = nanmean(junk,1);
    str = ['avg_pall.' xxname ' = xavg;'];
    eval(str);
  elseif xmm >  1 & xoo > 1
    xavg = squeeze(nanmean(junk,1));
    str = ['avg_pall.' xxname ' = xavg;'];
    eval(str);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix clouds as needed

if avg_pall.ctype > 001 & avg_pall.ctype < 101
  avg_pall.ctype = 101;
end
if avg_pall.ctype > 101 & avg_pall.ctype < 201
  avg_pall.ctype = 201;
end
if avg_pall.ctype > 201 & avg_pall.ctype < 301
  avg_pall.ctype = 301;
end

if avg_pall.ctype2 > 001 & avg_pall.ctype2 < 101
  avg_pall.ctype2 = 101;
end
if avg_pall.ctype2 > 101 & avg_pall.ctype2 < 201
  avg_pall.ctype2 = 201;
end
if avg_pall.ctype2 > 201 & avg_pall.ctype2 < 301
  avg_pall.ctype2 = 301;
end

avg_pall = fix_clouds_as_needed(avg_pall); %% see /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/

avg_pall.ctype(avg_pall.ctype >  100) = 201;
avg_pall.ctype(avg_pall.ctype <= 101) = 101;
avg_pall.ctype2(avg_pall.ctype2 >  100) = 201;
avg_pall.ctype2(avg_pall.ctype2 <= 101) = 101;

%%%%%%%%%%%%%%%%%%%%%%%%%

%% add in yy/mm/dd and salti/landfrac/emis

avg_pall.yy = floor(nanmean(avg_pall.yy)) * ones(size(avg_pall.stemp));
avg_pall.mm = floor(nanmean(avg_pall.mm)) * ones(size(avg_pall.stemp));
avg_pall.dd = floor(nanmean(avg_pall.dd)) * ones(size(avg_pall.stemp));
avg_pall.nlays = floor(avg_pall.nlays);
avg_pall.nlevs = floor(avg_pall.nlays) + 1;
%%% avg_pall.spres = a.pnew_op.spres;

sim_emis = load('/home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/tile_avg_salti_lf_emis.mat');

avg_pall.salti    = sim_emis.salti;
avg_pall.landfrac = sim_emis.landfrac;
avg_pall.nemis    = sim_emis.nemis;
avg_pall.rho      = sim_emis.rho;
avg_pall.emis     = sim_emis.emis;
avg_pall.efreq    = sim_emis.efreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now go through and check nlays vs spres
avg_pall.nlevs0 = avg_pall.nlevs;
avg_pall.nlays0 = avg_pall.nlays;
for ii = 1 : 4608
  junk = avg_pall.plevs(:,ii);
  junks = avg_pall.spres(ii);
  moo = find(junk <= junks);
  moo = max(moo);
  avg_pall.nlays(ii) = moo;
  avg_pall.nlevs(ii) = moo+1;
  avg_pall.blmult(ii) = (junk(moo)-junks)/(junk(moo)-junk(moo+1));
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(avg_pall.salti,72,64)') ,1), [-90 +90],[-180 +180]); title('salti'); colormap(jet)
figure(2); aslmap(2,rlat65,rlon73,smoothn((reshape(avg_pall.nemis,72,64)') ,1), [-90 +90],[-180 +180]); title('nemis'); colormap(jet)
figure(2); aslmap(2,rlat65,rlon73,smoothn((reshape(avg_pall.landfrac,72,64)') ,1), [-90 +90],[-180 +180]); title('landfrac'); colormap(jet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(avg_pall.stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 <ST>'); colormap(jet)
figure(2); aslmap(2,rlat65,rlon73,smoothn((reshape(avg_pall.mmw,72,64)') ,1),   [-90 +90],[-180 +180]); title('ERA5 <mmw>'); colormap(jet)
figure(3); junk = reshape(avg_pall.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(junk); colormap(jet); colorbar; caxis([200 300]); shading interp
figure(4); junk = reshape(avg_pall.gas_1,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(log10(junk)); colormap(jet); colorbar; caxis([15 21]); shading interp
figure(5); junk = reshape(avg_pall.gas_3,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(log10(junk)); colormap(jet); colorbar; caxis([14 18]); shading interp

figure(3); junk = reshape(avg_pall.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(meanvaluebin(latB2),log10(nanmean(avg_pall.plevs,2)),junk); colormap(jet); colorbar; caxis([200 300]); shading interp
figure(3); junk = reshape(avg_pall.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(meanvaluebin(latB2),nanmean(avg_pall.plevs,2),junk); colormap(jet); colorbar; caxis([200 300]); shading interp
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylim([10 1050])
figure(4); junk = reshape(avg_pall.gas_1,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(meanvaluebin(latB2),nanmean(avg_pall.plevs,2),log10(junk)); colormap(jet); colorbar; caxis([15 21]); shading interp
  set(gca,'yscale','linear'); set(gca,'ydir','reverse'); ylim([100 1050])  
figure(5); junk = reshape(avg_pall.gas_3,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(meanvaluebin(latB2),nanmean(avg_pall.plevs,2),log10(junk)); colormap(jet); colorbar; caxis([14 18]); shading interp
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylim([0.1 1050])  
