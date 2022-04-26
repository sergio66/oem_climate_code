% see ../FIND_NWP_MODEL_TRENDS/EvanManningCCRTrends/reader_12monthsloop.maddpath /home/sergio/MATLABCODE

%  The main field is the radiances, gridded at 10-degrees lat by 20-degrees lon but also with 10 angular bins and separation asc/desc and land/sea.

addpath /home/sergio/KCARTA/MATLAB

clear all; 

months = [1:12];
iCnt = 0;
for yy = 2002:2021
  mmS = 1; mmE = 12;
  if yy == 2002
    mmS = 09;
  elseif yy == 2021
    mmE = 08;
  end
  for mm = mmS : mmE
    iCnt = iCnt + 1;
    iaFound(iCnt) = 0;
    fname = ['/asl/s1/sergio/AIRS_L3/EVAN_MANNING_MONTHLY_CCR/stratrad_airs.cc.v7.' num2str(yy) '-' num2str(months(mm),'%02d') '.nc'];
    if exist(fname)
      iaFound(iCnt) = 1;
      saveyear(iCnt)  = yy;
      savemonth(iCnt) = months(mm);
      saveday(iCnt)   = 15;

      fprintf(1,'reading %3i %s \n',iCnt,fname);
      a = read_netcdf_lls(fname);
    
      %                           %% nchan satzen lon  lat L/O/A  Asc/Desc
      %size(a.l1b_airs.rad)       %% 2378     10   18   18   3       2
      %size(a.l1b_airs.rad_nobs)  %% 2378     10   18   18   3       2
    
      rads = squeeze(a.l1b_airs.rad(:,:,:,:,3,2));
      nobs  = squeeze(a.l1b_airs.rad_nobs(:,:,:,:,3,2));
    
      rads = squeeze(nanmean(rads,[2])); %% [2] = satzen
      nobs = squeeze(nansum(nobs,[2]));  %% [2] = satzen
    
      saverads(iCnt,:,:,:) = rads;
      savenobs(iCnt,:,:,:) = nobs;
    end
  end
end

error('opopo')

saverads = permute(saverads,[2 3 4 1]);
savenobs = permute(savenobs,[2 3 4 1]);

tCCR = rad2bt(a.l1b_airs.wnum,saverads);

tglobalCCR = squeeze(nanmean(tCCR,4));

pcolor(squeeze(tCCR(1291,:,:,1))'); colorbar; colormap jet
pcolor(squeeze(tglobalCCR(1291,:,:))'); colorbar; colormap jet

plot(squeeze(tCCR(1291,9,:,1))      %% chan x lon x lat x time
plot(squeeze(tglobalCCR(1291,9,:))) %% chan x lon x lat

iCntMax = iCnt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
days = (1:iCntMax)/12 * 001; %% so this is years
days = (1:iCntMax)/12 * 365; %% so this is days

yy = []; mm = []; dd = [];
for ii = 2002 : 2021
  clear yyx mmx ddx
  if ii == 2002
    inum = 4;
    yyx(1:inum) = ii;
    mmx = 9:12;
    ddx = ones(size(mmx)) * 15;
  elseif ii == 2021
    %inum = 7;
    %yyx(1:inum) = ii;
    %mmx = 1 : 7;
    inum = 8;
    yyx(1:inum) = ii;
    mmx = 1 : 8;
    ddx = ones(size(mmx)) * 15;
  else
    inum = 12;
    yyx(1:inum) = ii;
    mmx = 1:12;
    ddx = ones(size(mmx)) * 15;
  end
  fprintf(1,'%4i %2i \n',[ii inum])
  yy = [yy yyx];
  mm = [mm mmx];
  dd = [dd ddx];
end
rtime = utc2taiSergio(yy,mm,dd,ones(size(yy))*12.0);
dayOFtime = change2days(yy,mm,dd,2002);
days = dayOFtime;

addpath /home/sergio/MATLABCODE/FIND_TRENDS
warning off
for cc = 1 : 2378
  if mod(cc,1000) == 0
    fprintf(1,'+')
  elseif mod(cc,100) == 0
    fprintf(1,'.')
  end
  for ii = 1 : 18
    for jj = 1 : 18
      data = squeeze(double(tCCR(cc,ii,jj,:)));
      boo = find(isfinite(data));
      if length(boo) > 20
        [c d] = polyfit(days(boo)/365,data(boo),1);
        trend(cc,ii,jj) = c(1);
        const(cc,ii,jj) = c(2);
      else
        trend(cc,ii,jj) = NaN;
        const(cc,ii,jj) = NaN;
      end

      if length(boo) > 20
        [c d] = Math_tsfit_lin_robust(days(boo),data(boo),4);
        ccrtrend(cc,ii,jj) = c(2);
        ccrconst(cc,ii,jj) = c(1);
      else
        ccrtrend(cc,ii,jj) = NaN;
        ccrconst(cc,ii,jj) = NaN;
     end

    end
  end
end
warning on
fprintf(1,'\n');

%ccrtrend = bettertrend;
%ccrconst = betterconst;
comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_trend_ccr_evan_manning.m';
wnum2378 = a.l1b_airs.wnum;
save ccrtrends_2002_09_2021_08.mat comment trend const ccrtrend ccrconst wnum2378 comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% quick compare to our clearest quantile
compare_ccr_trends_against_Q16_L1c_trends

