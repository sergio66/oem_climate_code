%% close all; 
scr_siz = get(0,'ScreenSize');

figure(1); scr_siz = get(gcf); a0 = scr_siz.Position;  %% a0 = 1640         534         560         420

%% However, you forgot to specify the name of the ‘Position’ property. For example, here’s how to set the figure to be 500 pixels by 400 pixels:
%% set(gcf, 'Position',  [100+50, 100+50, 500, 400])

fac = 1.25;
figure(2); set(gcf, 'Position',  [100,    100, a0(3)*fac, a0(4)])
figure(3); set(gcf, 'Position',  [100+250, 100, a0(3)*fac, a0(4)])
figure(4); set(gcf, 'Position',  [100,    100+250, a0(3)*fac, a0(4)])
figure(5); set(gcf, 'Position',  [100+250, 100+250, a0(3)*fac, a0(4)])

