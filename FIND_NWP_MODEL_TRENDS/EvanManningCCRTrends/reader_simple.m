addpath /home/sergio/MATLABCODE

a=read_netcdf_lls('/asl/s1/sergio/AIRS_L3/EVAN_MANNING_MONTHLY_CCR/stratrad_airs.cc.v7.2003-09.nc');

a.l1b_airs.xtrack
                           %% nchan satzen lon  lat L/O/A  Asc/Desc
size(a.l1b_airs.rad)       %% 2378     10   18   18   3       2
size(a.l1b_airs.rad_nobs)  %% 2378     10   18   18   3       2

rads = squeeze(a.l1b_airs.rad(:,:,:,:,3,2));
nobs  = squeeze(a.l1b_airs.rad_nobs(:,:,:,:,3,2));

rads = squeeze(nanmean(rads,[2 3])); %% [2 3] = satzen/lon
nobs = squeeze(nansum(nobs,[2 3]));  %% [2 3] = satzen/lon

scatter(1:18,rad2bt(1231,rads(1291,:)),30,nobs(1291,:),'filled'); colormap jet; colorbar

sum(nobs(1291,:))
plot(a.l1b_airs.wnum,sum(nobs,2))

good = find(a.l1b_airs.nedn < 0.5);
plot(a.l1b_airs.wnum(good),sum(nobs(good,:),2))

allnobs  = squeeze(a.l1b_airs.rad_nobs(1291,:,:,:,3,:));
allnobs = nansum(allnobs(:))

i667 = find(a.l1b_airs.wnum >= 667,1)
allnobs  = squeeze(a.l1b_airs.rad_nobs(i667,:,:,:,3,:));
allnobs = nansum(allnobs(:))
