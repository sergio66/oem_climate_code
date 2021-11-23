goodtime = 1 : length(days); %% use ALL days
time  = days(goodtime);

[~,~,aauseQ,~] = size(save64x72_Q);
[~,~,aauseO3,~] = size(save64x72_O3);
[~,~,aauseT,~] = size(save64x72_T);

latbin = 20;
latbin = 32;
mon = 6;

inds = 1 : 12 : length(days);
inds = (inds-1) + mon;         %% once a year, not even smoothed ugh

for ii = 1 : 72
    stemp64x72_junk_anom(ii,:) = squeeze(save64x72_stemp(latbin,ii,inds)) - nanmean(squeeze(save64x72_stemp(latbin,ii,inds)));

    boo = squeeze(save64x72_Q(latbin,ii,:,inds));
    boo = nanmean(boo,2);
    boo = boo * ones(1,length(inds));
    woo = zeros(1,aauseQ,length(inds));
    woo(1,:,:) = boo;
    wah = squeeze(save64x72_Q(latbin,ii,:,inds));
    whos inds woo wah    
    Q64x72_junk_anom(ii,:,:)   = (squeeze(save64x72_Q(latbin,ii,:,inds)-woo)./woo);
end

    boo = squeeze(save64x72_O3(latbin,:,inds));
    boo = nanmean(boo,2);
    boo = boo * ones(1,length(inds));
    woo = zeros(1,aauseO3,length(inds));
    woo(1,:,:) = boo;
    whos inds woo    
    O364x72_junk_anom   = squeeze((save64x72_O3(latbin,:,inds)-woo)./woo);

    boo = squeeze(save64x72_T(latbin,:,inds));
    boo = nanmean(boo,2);
    boo = boo * ones(1,length(inds));
    woo = zeros(1,aauseT,length(inds));
    woo(1,:,:) = boo;
    whos inds woo    
    T64x72_junk_anom   = squeeze(save64x72_T(latbin,:,inds) - woo);

whos *junk_anom
figure(1); clf; plot(1:timespan,stemp64x72_junk_anom);
figure(2); clf; pcolor(1:timespan,1:12,Q64x72_junk_anom); shading interp; colormap jet; colorbar; title('Q anom');
figure(3); clf; pcolor(1:timespan,1:24,O364x72_junk_anom); shading interp; colormap jet; colorbar; title('WV anom');
figure(4); clf; pcolor(1:timespan,1:24,T64x72_junk_anom); shading interp; colormap jet; colorbar; title('T anom')
