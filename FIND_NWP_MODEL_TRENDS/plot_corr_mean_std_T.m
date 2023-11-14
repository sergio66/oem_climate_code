%% do 10 - 100, 100 - 400, 400 - 1000 mb

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10 - 100 mb
junk0 = miaow15(i10_100,:); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i10_100,:); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i10_100,:); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i10_100,:); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i10_100,:); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_global(:,1) = R(:,1); M_T_global(:,1) = m'; S_T_global(:,1) = s'; Slope_T_global(:,1) = linearfit(:,1); Frac_T_global(:,1) = frac_neg0pos(:,1);

junk0 = miaow15(i10_100,trend_rlat64_polar); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i10_100,trend_rlat64_polar); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i10_100,trend_rlat64_polar); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i10_100,trend_rlat64_polar); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i10_100,trend_rlat64_polar); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_polar(:,1) = R(:,1); M_T_polar(:,1) = m'; S_T_polar(:,1) = s'; Slope_T_polar(:,1) = linearfit(:,1); Frac_T_polar(:,1) = frac_neg0pos(:,1);

junk0 = miaow15(i10_100,trend_rlat64_midlat); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i10_100,trend_rlat64_midlat); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i10_100,trend_rlat64_midlat); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i10_100,trend_rlat64_midlat); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i10_100,trend_rlat64_midlat); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_midlat(:,1) = R(:,1); M_T_midlat(:,1) = m'; S_T_midlat(:,1) = s'; Slope_T_midlat(:,1) = linearfit(:,1); Frac_T_midlat(:,1) = frac_neg0pos(:,1);

junk0 = miaow15(i10_100,trend_rlat64_tropical); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i10_100,trend_rlat64_tropical); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i10_100,trend_rlat64_tropical); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i10_100,trend_rlat64_tropical); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i10_100,trend_rlat64_tropical); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_tropical(:,1) = R(:,1); M_T_tropical(:,1) = m'; S_T_tropical(:,1) = s'; Slope_T_tropical(:,1) = linearfit(:,1); Frac_T_tropical(:,1) = frac_neg0pos(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 100 - 400 mb
junk0 = miaow15(i100_400,:); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i100_400,:); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i100_400,:); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i100_400,:); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i100_400,:); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_global(:,2) = R(:,1); M_T_global(:,2) = m'; S_T_global(:,2) = s'; Slope_T_global(:,2) = linearfit(:,1); Frac_T_global(:,2) = frac_neg0pos(:,1);

junk0 = miaow15(i100_400,trend_rlat64_polar); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i100_400,trend_rlat64_polar); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i100_400,trend_rlat64_polar); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i100_400,trend_rlat64_polar); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i100_400,trend_rlat64_polar); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_polar(:,2) = R(:,1); M_T_polar(:,2) = m'; S_T_polar(:,2) = s'; Slope_T_polar(:,2) = linearfit(:,1); Frac_T_polar(:,2) = frac_neg0pos(:,1);

junk0 = miaow15(i100_400,trend_rlat64_midlat); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i100_400,trend_rlat64_midlat); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i100_400,trend_rlat64_midlat); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i100_400,trend_rlat64_midlat); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i100_400,trend_rlat64_midlat); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_midlat(:,2) = R(:,1); M_T_midlat(:,2) = m'; S_T_midlat(:,2) = s'; Slope_T_midlat(:,2) = linearfit(:,1); Frac_T_midlat(:,2) = frac_neg0pos(:,1);

junk0 = miaow15(i100_400,trend_rlat64_tropical); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i100_400,trend_rlat64_tropical); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i100_400,trend_rlat64_tropical); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i100_400,trend_rlat64_tropical); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i100_400,trend_rlat64_tropical); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_tropical(:,2) = R(:,1); M_T_tropical(:,2) = m'; S_T_tropical(:,2) = s'; Slope_T_tropical(:,2) = linearfit(:,1); Frac_T_tropical(:,2) = frac_neg0pos(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 400 - 1000 mb
junk0 = miaow15(i400_1000,:); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i400_1000,:); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i400_1000,:); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i400_1000,:); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i400_1000,:); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_global(:,3) = R(:,1); M_T_global(:,3) = m'; S_T_global(:,3) = s'; Slope_T_global(:,3) = linearfit(:,1); Frac_T_global(:,3) = frac_neg0pos(:,1);

junk0 = miaow15(i400_1000,trend_rlat64_polar); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i400_1000,trend_rlat64_polar); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i400_1000,trend_rlat64_polar); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i400_1000,trend_rlat64_polar); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i400_1000,trend_rlat64_polar); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_polar(:,3) = R(:,1); M_T_polar(:,3) = m'; S_T_polar(:,3) = s'; Slope_T_polar(:,3) = linearfit(:,1); Frac_T_polar(:,3) = frac_neg0pos(:,1);

junk0 = miaow15(i400_1000,trend_rlat64_midlat); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i400_1000,trend_rlat64_midlat); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i400_1000,trend_rlat64_midlat); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i400_1000,trend_rlat64_midlat); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i400_1000,trend_rlat64_midlat); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_midlat(:,3) = R(:,1); M_T_midlat(:,3) = m'; S_T_midlat(:,3) = s'; Slope_T_midlat(:,3) = linearfit(:,1); Frac_T_midlat(:,3) = frac_neg0pos(:,1);

junk0 = miaow15(i400_1000,trend_rlat64_tropical); junk0 = junk0(:);   j0 = nanmean(junk0); jj0 = nanstd(junk0);
  junk1 = miaow11(i400_1000,trend_rlat64_tropical); junk1 = junk1(:); j1 = nanmean(junk1); jj1 = nanstd(junk1);
  junk2 = miaow12(i400_1000,trend_rlat64_tropical); junk2 = junk2(:); j2 = nanmean(junk2); jj2 = nanstd(junk2);
  junk3 = miaow13(i400_1000,trend_rlat64_tropical); junk3 = junk3(:); j3 = nanmean(junk3); jj3 = nanstd(junk3);
  junk4 = miaow14(i400_1000,trend_rlat64_tropical); junk4 = junk4(:); j4 = nanmean(junk4); jj4 = nanstd(junk4);
  javg = nanmean([j0 j1 j2 j3 j4]); jjavg = nanmean([jj0 jj1 jj2 jj3 jj4]);
  boo = find(junk0 > javg - 2 * jjavg & junk0 < javg + 2 * jjavg & ...
             junk1 > javg - 2 * jjavg & junk1 < javg + 2 * jjavg & ...
             junk2 > javg - 2 * jjavg & junk2 < javg + 2 * jjavg & ...
             junk3 > javg - 2 * jjavg & junk3 < javg + 2 * jjavg & ...
             junk4 > javg - 2 * jjavg & junk4 < javg + 2 * jjavg);
  zall = [junk0(boo) junk1(boo) junk2(boo) junk3(boo) junk4(boo)]';
  if iBiasWRT_ERA5orUMBC > 0
    % iBiasWRT_ERA5orUMBC > = ==> wrt ERA5
    zall = zall([5 2 3 4 1],:);
  end
  wall = ones(size(junk0(boo)));
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
  R_T_tropical(:,3) = R(:,1); M_T_tropical(:,3) = m'; S_T_tropical(:,3) = s'; Slope_T_tropical(:,3) = linearfit(:,1); Frac_T_tropical(:,3) = frac_neg0pos(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%want a 1x3 tiled layout
ta = tiledlayout(3,5);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

dxyzaspect = [1.5 1.75 1]; %% not bad
dxyzaspect = [2 1.75 1]; %% not bad

smoothPts = 1;
smoothPts = 5;
smoothPts = 3;

iT = 0;

%% 10 - 100 mb
iT = iT + 1;
tafov(iT) = nexttile; imagesc([R_T_global(2:5,1) R_T_tropical(2:5,1) R_T_midlat(2:5,1) R_T_polar(2:5,1)]); colormap(tafov(iT),jet); caxis([0 +1]); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T CORRELATION'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [R_T_global(2:5,1) R_T_tropical(2:5,1) R_T_midlat(2:5,1) R_T_polar(2:5,1)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([Slope_T_global(2:5,1) Slope_T_tropical(2:5,1) Slope_T_midlat(2:5,1) Slope_T_polar(2:5,1)]); colormap(tafov(iT),usa2); caxis([-1 +1]); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T SLOPE'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [Slope_T_global(2:5,1) Slope_T_tropical(2:5,1) Slope_T_midlat(2:5,1) Slope_T_polar(2:5,1)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([M_T_global(2:5,1) M_T_tropical(2:5,1) M_T_midlat(2:5,1) M_T_polar(2:5,1)]); colormap(tafov(iT),usa2); caxis([-1 +1]*0.025); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T BIAS'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [M_T_global(2:5,1) M_T_tropical(2:5,1) M_T_midlat(2:5,1) M_T_polar(2:5,1)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([S_T_global(2:5,1) S_T_tropical(2:5,1) S_T_midlat(2:5,1) S_T_polar(2:5,1)]); colormap(tafov(iT),jet); caxis([0 +1]*0.025); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T STDDEV'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [S_T_global(2:5,1) S_T_tropical(2:5,1) S_T_midlat(2:5,1) S_T_polar(2:5,1)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([Frac_T_global(2:5,1) Frac_T_tropical(2:5,1) Frac_T_midlat(2:5,1) Frac_T_polar(2:5,1)]); colormap(tafov(iT),jet); caxis([0 +1]); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T +/-FRAC'); t.FontSize = 10; t.FontWeight = 'normal';
t = text(5,4,'10-100 mb','Rotation',90');
if iWriteCorrelNumbers > 0
  wah = [Frac_T_global(2:5,1) Frac_T_tropical(2:5,1) Frac_T_midlat(2:5,1) Frac_T_polar(2:5,1)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

%% 100 - 400 mb
iT = iT + 1;
tafov(iT) = nexttile; imagesc([R_T_global(2:5,2) R_T_tropical(2:5,2) R_T_midlat(2:5,2) R_T_polar(2:5,2)]); colormap(tafov(iT),jet); caxis([0 +1]); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  % t = title('T CORRELATION'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [R_T_global(2:5,2) R_T_tropical(2:5,2) R_T_midlat(2:5,2) R_T_polar(2:5,2)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([Slope_T_global(2:5,2) Slope_T_tropical(2:5,2) Slope_T_midlat(2:5,2) Slope_T_polar(2:5,2)]); colormap(tafov(iT),usa2); caxis([-1 +1]); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  % t = title('T SLOPE'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [Slope_T_global(2:5,2) Slope_T_tropical(2:5,2) Slope_T_midlat(2:5,2) Slope_T_polar(2:5,2)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([M_T_global(2:5,2) M_T_tropical(2:5,2) M_T_midlat(2:5,2) M_T_polar(2:5,2)]); colormap(tafov(iT),usa2); caxis([-1 +1]*0.025); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  % t = title('T BIAS'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [M_T_global(2:5,2) M_T_tropical(2:5,2) M_T_midlat(2:5,2) M_T_polar(2:5,2)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([S_T_global(2:5,2) S_T_tropical(2:5,2) S_T_midlat(2:5,2) S_T_polar(2:5,2)]); colormap(tafov(iT),jet); caxis([0 +1]*0.025); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  % t = title('T STDDEV'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [S_T_global(2:5,2) S_T_tropical(2:5,2) S_T_midlat(2:5,2) S_T_polar(2:5,2)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([Frac_T_global(2:5,2) Frac_T_tropical(2:5,2) Frac_T_midlat(2:5,2) Frac_T_polar(2:5,2)]); colormap(tafov(iT),jet); caxis([0 +1]); 
  %colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %% t = title('T +/-FRAC'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [Frac_T_global(2:5,2) Frac_T_tropical(2:5,2) Frac_T_midlat(2:5,2) Frac_T_polar(2:5,2)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

t = text(5,4,'100-400 mb','Rotation',90');

%% 400 - 1000 mb
iT = iT + 1;
tafov(iT) = nexttile; imagesc([R_T_global(2:5,3) R_T_tropical(2:5,3) R_T_midlat(2:5,3) R_T_polar(2:5,3)]); colormap(tafov(iT),jet); caxis([0 +1]); 
  colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %% t = title('T CORRELATION'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [R_T_global(2:5,3) R_T_tropical(2:5,3) R_T_midlat(2:5,3) R_T_polar(2:5,3)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([Slope_T_global(2:5,3) Slope_T_tropical(2:5,3) Slope_T_midlat(2:5,3) Slope_T_polar(2:5,3)]); colormap(tafov(iT),usa2); caxis([-1 +1]); 
  colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %% t = title('T SLOPE'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [Slope_T_global(2:5,3) Slope_T_tropical(2:5,3) Slope_T_midlat(2:5,3) Slope_T_polar(2:5,3)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([M_T_global(2:5,3) M_T_tropical(2:5,3) M_T_midlat(2:5,3) M_T_polar(2:5,3)]); colormap(tafov(iT),usa2); caxis([-1 +1]*0.025); 
  colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %% t = title('T BIAS'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [M_T_global(2:5,3) M_T_tropical(2:5,3) M_T_midlat(2:5,3) M_T_polar(2:5,3)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([S_T_global(2:5,3) S_T_tropical(2:5,3) S_T_midlat(2:5,3) S_T_polar(2:5,3)]); colormap(tafov(iT),jet); caxis([0 +1]*0.025); 
  colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %% t = title('T STDDEV'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [S_T_global(2:5,3) S_T_tropical(2:5,3) S_T_midlat(2:5,3) S_T_polar(2:5,3)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end

iT = iT + 1;
tafov(iT) = nexttile; imagesc([Frac_T_global(2:5,3) Frac_T_tropical(2:5,3) Frac_T_midlat(2:5,3) Frac_T_polar(2:5,3)]); colormap(tafov(iT),jet); caxis([0 +1]); 
  colorbar('southoutside');
  set(gca,'ytick',[1:4],'yticklabel',modelnames(2:5))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'}); 
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %% t = title('T +/-FRAC'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = [Frac_T_global(2:5,3) Frac_T_tropical(2:5,3) Frac_T_midlat(2:5,3) Frac_T_polar(2:5,3)]';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = wah(iii,jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.10,jjj+0.25,junkstr,'Fontsize',7,'Rotation',90);
    end
  end
end
t = text(5,4,'400-1000 mb','Rotation',90');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 3 4 5   7 8 9 10   12 13 14 15]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end

% Remove all xtick labels except for 2nd row
for ii = [1 2 3 4 5   6 7 8 9 10]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end
