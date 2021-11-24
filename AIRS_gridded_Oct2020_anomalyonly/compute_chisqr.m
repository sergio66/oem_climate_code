for jj = 15 : 26
  jj
  for dd = 1 : 365
    fname = ['SAVE_LW_noCFC11_noN2O/OutputAnomaly_OBS/' num2str(jj) '/anomtest_timestep' num2str(dd) '.mat'];
    a = load(fname);
    obs = a.rateset.rates;
    unc = a.rateset.unc_rates;
    cal = a.oem.fit';
    chans = a.jacobian.chanset;
    nonzero = find(isfinite(obs) & isfinite(cal) & isfinite(unc) & unc > eps);
    boo = (obs-cal)./unc;
    allchisqr(jj-15+1,dd) = (boo(nonzero)'*boo(nonzero))/length(nonzero);
    chanschisqr(jj-15+1,dd) = (boo(chans)'*boo(chans))/length(chans);
    %if mod(dd,100) == 0
    %  figure(1); plot(allchisqr,'.-');
    %   figure(2); plot(chanschisqr,'.-');
    %  pause(0.1)
    %end
  end
end

figure(1); plot(allchisqr(:),'.-');
figure(2); plot(chanschisqr(:),'.-');
figure(3); plot(0:1000:100000,histc(allchisqr(:),0:1000:100000),'b',0:10:1000,histc(chanschisqr(:),0:10:1000),'r')
xlim([0 20000])
figure(3); plot(0:10:1000,histc(chanschisqr(:),0:10:1000),'r',0:10:1000,smooth(histc(chanschisqr(:),0:10:1000),10),'m')
