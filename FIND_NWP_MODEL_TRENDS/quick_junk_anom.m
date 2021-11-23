goodtime = 1 : length(days); %% use ALL days
time  = days(goodtime);

[x,aauseQ,x] = size(save_Q);
[x,aauseO3,x] = size(save_O3);
[x,aauseT,x] = size(save_T);


latbin = 20;
mon = 6;

    inds = 1 : 12 : length(days);
    inds = (inds-1) + mon;
    stemp_junk_anom = save_stemp(latbin,inds) - nanmean(save_stemp(latbin,inds));

    boo = squeeze(save_Q(latbin,:,inds));
    boo = nanmean(boo,2);
    boo = boo * ones(1,length(inds));
    woo = zeros(1,aauseQ,length(inds));
    woo(1,:,:) = boo;
    wah = save_Q(latbin,:,inds);
    whos inds woo wah    
    Q_junk_anom   = squeeze((save_Q(latbin,:,inds)-woo)./woo);

    boo = squeeze(save_O3(latbin,:,inds));
    boo = nanmean(boo,2);
    boo = boo * ones(1,length(inds));
    woo = zeros(1,aauseO3,length(inds));
    woo(1,:,:) = boo;
    whos inds woo    
    O3_junk_anom   = squeeze((save_O3(latbin,:,inds)-woo)./woo);

    boo = squeeze(save_T(latbin,:,inds));
    boo = nanmean(boo,2);
    boo = boo * ones(1,length(inds));
    woo = zeros(1,aauseT,length(inds));
    woo(1,:,:) = boo;
    whos inds woo    
    T_junk_anom   = squeeze(save_T(latbin,:,inds) - woo);

whos *junk_anom
figure(1); clf; plot(1:timespan,stemp_junk_anom);
figure(2); clf; pcolor(1:timespan,1:12,Q_junk_anom); shading interp; colormap jet; colorbar; title('Q anom');
figure(3); clf; pcolor(1:timespan,1:24,O3_junk_anom); shading interp; colormap jet; colorbar; title('WV anom');
figure(4); clf; pcolor(1:timespan,1:24,T_junk_anom); shading interp; colormap jet; colorbar; title('T anom')
