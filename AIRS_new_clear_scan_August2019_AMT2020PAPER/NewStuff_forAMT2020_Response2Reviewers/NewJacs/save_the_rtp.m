addpath /asl/matlib/h4tools

%{
pbefore.plevs = nanmean(poem.plevs,2) * ones(1,length(pbefore.stemp));
pafter.plevs = nanmean(poem.plevs,2) * ones(1,length(pafter.stemp));

pbefore.palts = nanmean(poem.palts,2) * ones(1,length(pbefore.stemp));
pafter.palts = nanmean(poem.palts,2) * ones(1,length(pafter.stemp));

pbefore.zobs    = 705000 * ones(size(pbefore.stemp));    pafter.zobs    = 705000 * ones(size(pafter.stemp));
pbefore.salti   = 0.0 * ones(size(pbefore.stemp));       pafter.salti   = 0.0 * ones(size(pafter.stemp));
pbefore.scanang = 22.0 * ones(size(pbefore.stemp));      pafter.scanang = 22.0 * ones(size(pafter.stemp));
pbefore.satzen = 22.0 * ones(size(pbefore.stemp));       pafter.satzen = 22.0 * ones(size(pafter.stemp));
pbefore.nlevs = 98 * ones(size(pbefore.stemp));        pafter.nlevs = 98 * ones(size(pafter.stemp));

pbefore.nemis = 19 * ones(size(pbefore.stemp));        pafter.nemis = 19 * ones(size(pafter.stemp));
pbefore.efreq = nanmean(poem.efreq,2) * ones(1,length(pbefore.stemp));
pbefore.emis  = nanmean(poem.emis,2) * ones(1,length(pbefore.stemp));
pbefore.rho   = (1-pbefore.emis)/pi;
pafter.efreq  = nanmean(poem.efreq,2) * ones(1,length(pafter.stemp));
pafter.emis   = nanmean(poem.emis,2) * ones(1,length(pafter.stemp));
pafter.rho    = (1-pafter.emis)/pi;
%}

%rtpwrite('pbefore.op.rtp',hx,[],pbefore,[]);  %% done Mar 30, 2020
%rtpwrite('pafter.op.rtp',hx,[],pafter,[]);    %% done Mar 30, 2020
rtpwrite('pbefore2.op.rtp',hx,[],pbefore,[]);  %% done Apr 23, 2020
rtpwrite('pafter2.op.rtp',hx,[],pafter,[]);    %% done Apr 23, 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xnameB = fieldnames(pbefore);
xnameA = fieldnames(pafter);
for ll = 1 : 40
  for ii = 1 : length(xnameB)
    boo = find(pbefore.latbin == ll);
    str = ['junk = pbefore.' xnameB{ii} ';'];
    eval(str)
    [mm,nn] = size(junk);
    if mm == 1
      str = ['pbeforeavg.' xnameB{ii} '(ll) = nanmean(junk(boo));'];
    else
      str = ['pbeforeavg.' xnameB{ii} '(:,ll) = nanmean(junk(:,boo),2);'];
    end
    eval(str)
  end

  for ii = 1 : length(xnameA)
    boo = find(pafter.latbin == ll);
    str = ['junk = pafter.' xnameA{ii} ';'];
    eval(str)
    [mm,nn] = size(junk);
    if mm == 1
      str = ['pafteravg.' xnameA{ii} '(ll) = nanmean(junk(boo));'];
    else
      str = ['pafteravg.' xnameA{ii} '(:,ll) = nanmean(junk(:,boo),2);'];
    end
    eval(str)
  end
end

figure(1); plot(pbefore.latbin,pbefore.stemp,'o-',pbeforeavg.latbin,pbeforeavg.stemp,'ro-')
figure(2); pcolor(pbeforeavg.ptemp); colorbar; caxis([200 300]); colormap jet
  set(gca,'ydir','reverse'); title('avg ptemp')
figure(3); pcolor(log10(pbeforeavg.gas_1)); colorbar; caxis(log10([1e15 1e22])); colormap jet
  set(gca,'ydir','reverse'); title('avg WV')
figure(4); pcolor(log10(pbeforeavg.gas_3)); colorbar; caxis(log10([1e12 1e18])); colormap jet
  set(gca,'ydir','reverse'); title('avg O3')

%rtpwrite('pbeforeavg.op.rtp',hx,[],pbeforeavg,[]);  %% done Mar 30, 2020
%rtpwrite('pafteravg.op.rtp',hx,[],pafteravg,[]);    %% done Mar 30, 2020
rtpwrite('pbeforeavg2.op.rtp',hx,[],pbeforeavg,[]);  %% done Apr 23, 2020
rtpwrite('pafteravg2.op.rtp',hx,[],pafteravg,[]);    %% done Apr 23, 2020

figure(1); pcolor(pafteravg.diff_mean_ptemp); shading interp; colorbar; set(gca,'ydir','reverse'); title('mean(Tz after-Tz before)');
figure(2); pcolor(pafteravg.diff_std_ptemp); shading interp; colorbar; set(gca,'ydir','reverse'); title('std(Tz after-Tz before)');
figure(3); pcolor(pafteravg.diff_mean_gas_1); shading interp; colorbar; set(gca,'ydir','reverse'); title('mean(WV after./WV before)');
figure(4); pcolor(pafteravg.diff_std_gas_1); shading interp; colorbar; set(gca,'ydir','reverse'); title('std(WV after./WV before)');
figure(5); pcolor(pafteravg.diff_mean_gas_3); shading interp; colorbar; set(gca,'ydir','reverse'); title('mean(O3 after./O3 before)');
figure(6); pcolor(pafteravg.diff_std_gas_3); shading interp; colorbar; set(gca,'ydir','reverse'); title('std(O3 after./O3 before)');

figure(1); caxis([-0.01 +0.01])
figure(2); caxis([0 +0.05])
