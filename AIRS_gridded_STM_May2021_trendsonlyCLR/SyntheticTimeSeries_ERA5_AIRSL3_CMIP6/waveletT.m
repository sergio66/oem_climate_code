% https://www.mathworks.com/help/wavelet/time-frequency-analysis.html
% https://www.mathworks.com/help/wavelet/ug/practical-introduction-to-continuous-analysis.html

%% junk = A sin (3 2 pi t/365) = A sin (6 pi t') = A sin (w t'); thus T = 2 * pi/w = 2 pi / 6 pi = 1/3 == > f = 1/T = 3
junk = 2*cos(1*2*pi*dayOFtime/365) + 3*sin(3*2*pi*dayOFtime/365); junk = junk + 0.1*randn(size(junk)); plot(junk)
dt = 1/12; fs = 1/dt; %% once a month, so frequency sampling = 12 times per year;

[cfs,fff] = cwt(junk,fs);
%helperHyperbolicChirpPlot(cfs,fff,2002+dayOFtime/365);
imagesc(2002+dayOFtime/365,fff,abs(cfs))
xlabel('Time')
ylabel('Frequency (yr-1)')
axis xy
caxis([0.025 0.25])
title('CWT of ECG Data')

pspectrum(junk,19*12,'spectrogram'); colormap jet;

wcoherence(stempjunk,tcalc(1520,:),12,'PhaseDisplayThreshold',0.7)
title('Wavelet coherence SKT and BT1231'); ylabel('Freq 1/yr'); xlabel('Time yr');

figure(2); plot(dayOFtime/365,stempjunk-mean(stempjunk),dayOFtime/365,tcalc(1520,:)-mean(tcalc(1520,:)))

figure(3);
wcoherence(stempjunk-mean(stempjunk),tcalc(1520,:)-mean(tcalc(1520,:),years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(stempjunk-mean(stempjunk),tcalc(1520,:)-mean(tcalc(1520,:)),years(1/12),'PhaseDisplayThreshold',0.7);
title('Wavelet coherence SKT and BT1231'); ylabel('Time Period yr'); xlabel('Time yr');

