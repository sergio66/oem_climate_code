% https://www.mathworks.com/help/signal/ref/hht.html
%% junk = A sin (3 2 pi t/365) = A sin (6 pi t') = A sin (w t'); thus T = 2 * pi/w = 2 pi / 6 pi = 1/3 == > f = 1/T = 3
junk = 2*cos(1*2*pi*dayOFtime/365) + 3*sin(3*2*pi*dayOFtime/365); junk = junk + 0.1*randn(size(junk)); plot(junk)
dt = 1/12; fs = 1/dt; %% once a month, so frequency sampling = 12 times per year;
emd(junk);
imf = emd(junk,'Display',1);
figure(2); [x,f,t,imfinsf,imfinse] = hht(imf,fs); hht(imf,fs); colormap(jet); xlabel('Time (years)'); ylabel('Frequency yr-1');
dfs = fs/2/100; fsgrid = (0:100)*dfs; % https://www.mathworks.com/help/signal/ref/hht.html
figure(3); pcolor(dayOFtime/365+2002,fsgrid,log10(full(x))); colormap(jet); colorbar; shading flat; %% x is sparse, so convert to full
plot(t,imfinsf)

%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/12; fs = 1/dt; %% once a month, so frequency sampling = 12 times per year;
emd(stempjunk);
imf = emd(stempjunk,'Display',1);
figure(5); [x,f,t,imfinsf,imfinse] = hht(imf,fs); colormap(jet); xlabel('SURFACE TEMP Time (years)'); ylabel('Frequency yr-1');
dfs = fs/2/100; fsgrid = (0:100)*dfs; % https://www.mathworks.com/help/signal/ref/hht.html
figure(6); pcolor(dayOFtime/365+2002,fsgrid,full(x)); colormap(jet); colorbar; shading flat; %% x is sparse, so convert to full
b = Math_tsfit_lin_robust(dayOFtime,stempjunk,4); b(2)

dt = 1/12; fs = 1/dt; %% once a month, so frequency sampling = 12 times per year;
emd(tcalc(1520,:));
imf = emd(tcalc(1520,:),'Display',1);
figure(8); [x,f,t,imfinsf,imfinse] = hht(imf,fs); colormap(jet); xlabel('BT 1231 Time (years)'); ylabel('Frequency yr-1');
dfs = fs/2/100; fsgrid = (0:100)*dfs; % https://www.mathworks.com/help/signal/ref/hht.html
figure(9); pcolor(dayOFtime/365+2002,fsgrid,full(x)); colormap(jet); colorbar; shading flat; %% x is sparse, so convert to full
b = Math_tsfit_lin_robust(dayOFtime,tcalc(1520,:),4); b(2)
