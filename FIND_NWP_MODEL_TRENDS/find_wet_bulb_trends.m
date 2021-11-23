%% see /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/wet_bulb_layer_dew_point_temperature
CtoK = 273.13;
Tsurf  = p.stemp - CtoK;  %% in Centigrade
RHsurf = RH0(97,:);  %%%% ORIG
RHsurf = RHSurf0;    %%%% HA, BETTER
Tw = Tsurf .* atan(0.151977 .* (RHsurf + 8.313659).^(1/2)) + atan(Tsurf + RHsurf) - atan(RHsurf - 1.676331) + 0.00391838*(RHsurf).^(3/2) .* atan(0.023101 * RHsurf) - 4.686035; %% in Centograde
Tw = Tw + CtoK;  %% in K

Tsurfpert  = pert.stemp - CtoK;  %% in Centigrade
RHsurfpert = RHpert(97,:);  %%%% ORIG
RHsurfpert = RHSurfpert;    %%%% HA, BETTER
Twpert = Tsurfpert .* atan(0.151977 .* (RHsurfpert + 8.313659).^(1/2)) + atan(Tsurfpert + RHsurfpert) - atan(RHsurfpert - 1.676331) + 0.00391838*(RHsurfpert).^(3/2) .* atan(0.023101 * RHsurfpert) - 4.686035; %% in Centograde
Twpert = Twpert + CtoK;  %% in K

dRHsurf = [RHSurfpert-RHSurf0];
figure(33); aslmap(33,rlat65,rlon73,maskLFmatr.*smoothn((reshape(RHsurf,72,64)'),1),[-90 +90],[-180 +180]);   colormap(jet);
figure(33); aslmap(33,rlat65,rlon73,maskLFmatr.*smoothn((reshape(dRHsurf,72,64)'),1),[-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.5 +0.5])

%dTw/dt = dTw/dTs dTs/dt + dTw/dRH dRH/dt     d/dx(atan(x)) = 1/(1+x^2)
dTw_dTs = atan(0.151977 .* (RHsurf + 8.313659).^(1/2)) + 1./(1+(Tsurf + RHsurf).^2);
ajunk = 1./((0.151977 .* (RHsurf + 8.313659).^(1/2)).^2 + 1) * 0.151977/2./sqrt(RHsurf + 8.313659); ajunk = Tsurf .* (ajunk); 
bjunk = (Tsurf + RHsurf).^2;     bjunk = 1./(1+bjunk);
cjunk = (RHsurf - 1.676331).^2;  cjunk = 1./(1+cjunk);
djunk = 0.00391838 * (3/2*(RHsurf).^(1/2).* atan(0.023101 * RHsurf) + (RHsurf).^(3/2) ./ (1 + (0.023101 * RHsurf).^2) * 0.023101);
dTw_dRH = ajunk + bjunk - cjunk + djunk;
dTw = dTw_dTs .* results(:,6)' + dTw_dRH .* deltaRH(97,:); %%%% ORIG
dTw = dTw_dTs .* results(:,6)' + dTw_dRH .* dRHsurf;       %%%% BETTER
dTw = Twpert - Tw;                                         %%%% HA, EVEN BETTER

figure(33); aslmap(33,rlat65,rlon73,maskLFmatr.*smoothn((reshape(p.stemp,72,64)')-CtoK,1),[-90 +90],[-180 +180]);   colormap(jet);  title('ST');        caxis([220 320]-CtoK); caxis([-40 +40])
figure(34); aslmap(34,rlat65,rlon73,maskLFmatr.*smoothn((reshape(RH0(97,:),72,64)'),1),[-90 +90],[-180 +180]);      colormap(jet);  title('surf RH0');  caxis([50 100])

figure(35); aslmap(35,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,6),72,64)'),1),[-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt ST');        caxis([-0.1 +0.1])
%figure(36); aslmap(36,rlat65,rlon73,maskLFmatr.*smoothn((reshape(deltaRH(97,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt surf RH');  caxis([-0.5 +0.5])
figure(36); aslmap(36,rlat65,rlon73,maskLFmatr.*smoothn((reshape(RHSurfpert-RHSurf0,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt surf RH');  caxis([-0.5 +0.5])
%junk = smoothn((reshape(deltaRH(97,:),72,64)'),1); junk = sign(junk).*log10(abs(junk));  aslmap(36,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]); colormap(llsmap5);  title('log (d/dt surf RH) ');  caxis([-0.5 +0.5])

figure(37); aslmap(37,rlat65,rlon73,maskLFmatr.*smoothn((reshape(Tw,72,64)')-CtoK,1),[-90 +90],[-180 +180]);        colormap(jet);  title('ST_{wet}');  caxis([220 320]-CtoK); caxis([-40 +40])
figure(38); aslmap(38,rlat65,rlon73,maskLFmatr.*smoothn((reshape(dTw,72,64)'),1),[-90 +90],[-180 +180]);            colormap(llsmap5);  title('d/dt ST_{wet}');  caxis([-0.1 +0.1])

%% now find time to get to 35 C
TwC = Tw - CtoK;
time2apocalypse = (35 - TwC); 
figure(39); aslmap(39,rlat65,rlon73,maskLFmatr.*smoothn((reshape(time2apocalypse,72,64)'),1),[-90 +90],[-180 +180]);            colormap(jet);  title('ST_{wet}-35C');  caxis([-10 50])

time2apocalypse = (35 - TwC)./dTw; 
time2apocalypse(time2apocalypse < 0) = NaN;
figure(40); aslmap(40,rlat65,rlon73,maskLFmatr.*smoothn((reshape(time2apocalypse,72,64)'),1),[-90 +90],[-180 +180]);            colormap(jet);  title('Years to Apocalype');  caxis([-1 1000])
figure(40); aslmap(40,rlat65,rlon73,maskLFmatr.*smoothn((reshape(log10(time2apocalypse),72,64)'),1),[-90 +90],[-180 +180]);     colormap(jet);  title('log10(Years to Apocalype)');  caxis([0 5])

waboo = time2apocalypse;                               wamoo = find(isnan(time2apocalypse));
bonk = reshape(time2apocalypse,72,64)';                wonk = find(isnan(time2apocalypse));
sonk = smoothn((reshape(time2apocalypse,72,64)'),1);   sonk(sonk < 0) = NaN; sonk(wonk) = NaN;
figure(40); aslmap(40,rlat65,rlon73,maskLFmatr.*sonk,[-90 +90],[-180 +180]);            colormap(jet);  title('Years to Apocalype');         caxis([-1 1000])
figure(40); aslmap(40,rlat65,rlon73,maskLFmatr.*log10(sonk),[-90 +90],[-180 +180]);     colormap(jet);  title('log10(Years to Apocalype)');  caxis([-5 10]); caxis([0 5]); colormap(hot)
