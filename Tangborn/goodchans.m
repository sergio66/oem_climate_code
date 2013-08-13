function  goodchannels1 = goodchan(f,robs,goodchannels);

clf
fg = f(goodchannels);
tg = rad2bt(f(goodchannels),robs(goodchannels));

tdiff = diff(tg);
ii = find(abs(tdiff) < 3);
goodchannels1 = goodchannels(ii);

fg = f(goodchannels1);
tg = rad2bt(f(goodchannels1),robs(goodchannels1));

subplot(211); plot(fg,tg)

mn = mean(diff(abs(tg)))
subplot(212); plot(fg(1:length(fg)-1),diff(abs(tg)))