addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB

iCnt = 0;
for yy = 2003:2021
  iCnt = iCnt + 1;
  a=read_netcdf_lls(['/asl/s1/sergio/AIRS_L3/EVAN_MANNING_MONTHLY_CCR/stratrad_airs.cc.v7.' num2str(yy) '-07.nc']);

  %                           %% nchan satzen lon  lat L/O/A  Asc/Desc
  %size(a.l1b_airs.rad)       %% 2378     10   18   18   3       2
  %size(a.l1b_airs.rad_nobs)  %% 2378     10   18   18   3       2

  rads = squeeze(a.l1b_airs.rad(:,:,:,:,3,2));
  nobs  = squeeze(a.l1b_airs.rad_nobs(:,:,:,:,3,2));

  rads = squeeze(nanmean(rads,[2 3])); %% [2 3] = satzen/lon
  nobs = squeeze(nansum(nobs,[2 3]));  %% [2 3] = satzen/lon

  saverads(iCnt,:,:) = rads;
  savenobs(iCnt,:,:) = nobs;
end

saverads = permute(saverads,[2 3 1]);
savenobs = permute(savenobs,[2 3 1]);

tCCR = rad2bt(a.l1b_airs.wnum,saverads);

tglobalCCR = squeeze(nanmean(tCCR,2));
for ii = 1 : 2378
  [c d]= polyfit(1:19,tglobalCCR(ii,:),1);
  trend(ii) = c(1);
  const(ii) = c(2);
end
plot(a.l1b_airs.wnum,const); title('July mean BTobs');
plot(a.l1b_airs.wnum,trend); title('July trend mean BTobs');

plotaxis2; axis([640 1640 -0.1 +0.05])
