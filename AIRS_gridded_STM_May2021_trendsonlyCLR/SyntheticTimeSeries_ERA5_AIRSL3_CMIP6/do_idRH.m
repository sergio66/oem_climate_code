%  idRH = +1;   %% keep WV   constant, no clouds
%  idRH = +2;   %% keep CO2  constant, no clouds
%  idRH = +3;   %% keep T,ST constant, no clouds
%  idRH = +4;   %% keep RH   constant, no clouds
%  idRH = +5;   %% put in everything, including clouds
%  idRH = +6;   %% increase RH by 0.01, no clouds
%  idRH = +7;   %% do nothing, no clouds

if idRH == 1
  ind = (1:72);
  pconst  = p72.plevs(:,ind);
  stconst = p72.stemp(ind);
  tconst = p72.ptemp(:,ind);
  wvconst = p72.gas_1(:,ind);
  rhconst = convert_humidity (pconst,tconst,wvconst,'specific humidity','relative humidity');

  wvconst = p72.gas_1(:,ind);
  for iTime = 1 : length(p72.stemp)/72
    indTime = (iTime-1)*72 + (1:72);
    p72.gas_1(:,indTime) = wvconst;
    rhnew = convert_humidity (p72.plevs(:,indTime),p72.ptemp(:,indTime),p72.gas_1(:,indTime),'specific humidity','relative humidity');
    drhsave(iTime) = rhnew(20,36)-rhconst(20,36);
  end
  junk = p72.gas_1(:,ind+72);
  figure(1); clf; pcolor(wvconst-junk); shading interp; colormap(usa2); colorbar; title('idRH = 1 ... wv const'); set(gca,'ydir','reverse')
  figure(6); clf; plot(drhsave); title(['idRH = ' num2str(idRH) ' \delta RH at center of 37 x 72 set']); set(gca,'fontsize',10)
    pause(0.1)

elseif idRH == 2
  ind = (1:72);
  pconst  = p72.plevs(:,ind);
  stconst = p72.stemp(ind);
  tconst = p72.ptemp(:,ind);
  wvconst = p72.gas_1(:,ind);
  rhconst = convert_humidity (pconst,tconst,wvconst,'specific humidity','relative humidity');

  co2const = p72.co2ppm(ind);
  for iTime = 1 : length(p72.stemp)/72
    indTime = (iTime-1)*72 + (1:72);
    p72.co2ppm(indTime) = co2const;
    rhnew = convert_humidity (p72.plevs(:,indTime),p72.ptemp(:,indTime),p72.gas_1(:,indTime),'specific humidity','relative humidity');
    drhsave(iTime) = rhnew(20,36)-rhconst(20,36);
  end
  figure(1); clf; plot(p72.co2ppm); title('idRH = 2 ... co2 const')
  figure(6); clf; plot(drhsave); title(['idRH = ' num2str(idRH) ' \delta RH at center of 37 x 72 set']); set(gca,'fontsize',10)
    pause(0.1)

elseif idRH == 3
  ind = (1:72);
  pconst  = p72.plevs(:,ind);
  stconst = p72.stemp(ind);
  tconst = p72.ptemp(:,ind);
  wvconst = p72.gas_1(:,ind);
  rhconst = convert_humidity (pconst,tconst,wvconst,'specific humidity','relative humidity');

  stconst = p72.stemp(ind);
  tconst = p72.ptemp(:,ind);
  for iTime = 1 : length(p72.stemp)/72
    indTime = (iTime-1)*72 + (1:72);
    p72.stemp(indTime) = stconst;
    p72.ptemp(:,indTime) = tconst;
    rhnew = convert_humidity (p72.plevs(:,indTime),p72.ptemp(:,indTime),p72.gas_1(:,indTime),'specific humidity','relative humidity');
    drhsave(iTime) = rhnew(20,36)-rhconst(20,36);
  end
  junk = p72.ptemp(:,ind+72);
  figure(1); clf; pcolor(tconst-junk); shading interp; colormap(usa2); colorbar; title('idRH = 2 ... T const'); set(gca,'ydir','reverse')
  figure(6); clf; plot(drhsave); title(['idRH = ' num2str(idRH) ' \delta RH at center of 37 x 72 set']); set(gca,'fontsize',10)
    pause(0.1)

elseif idRH == 4 | idRH == 6
  %% klayers
  %         20   mass mixing ratio in (g/kg), dry air
  %              Grams of gas X per kilogram of "dry air"
  %
  %         21   mass mixing ratio in (g/g) or (kg/kg), dry air
  %              Grams of gas X per gram of "dry air"
  %
  % while "convert_humidity" requires mixing ratio in kg/kg (not g/kg)
  initialWV = p72.gas_1;
  ind = (1:72);
  pconst  = p72.plevs(:,ind);
  stconst = p72.stemp(ind);
  tconst = p72.ptemp(:,ind);
  wvconst = p72.gas_1(:,ind);
  rhconst = convert_humidity (pconst,tconst,wvconst,'specific humidity','relative humidity');

  %% look at trend paper : for const RH : Eqn 4 : de = e (Lv/Rv dT/T^2)    where Lv = 
  %% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/guess_wv_surface.. for Lo,Rv
  Lo = 2.5e6;  %%% J/kg
  Rv = 461.52; %%% J/kg/K
  for iTime = 1 : length(p72.stemp)/72
    indTime = (iTime-1)*72 + (1:72);
    p72.stemp(indTime) = p72.stemp(indTime);
    p72.ptemp(:,indTime) = p72.ptemp(:,indTime);
    dT = p72.ptemp(:,indTime) - tconst;
    q  = p72.gas_1(:,indTime);
    dq = wvconst .* (Lo/Rv * dT ./p72.ptemp(:,indTime) ./p72.ptemp(:,indTime));
    if idRH == 4 & iTime > 1
      %% keep RH constant from time t = 0
      p72.gas_1(:,indTime) = wvconst + dq;
    elseif idRH == 6 & iTime > 1
      %% increase RH time t = 0
      RHincrease = 0.01; %% 1 perecent increase per year, divide by 12 per month

      %% see https://vortex.plymouth.edu/~stmiller/stmiller_content/Publications/AtmosRH_Equations_Rev.pdf4
      rh = convert_humidity (p72.plevs(:,indTime),p72.ptemp(:,indTime),p72.gas_1(:,indTime),'specific humidity','relative humidity');

      %% esat(373) = 1001 mb, yay
      satvp = esat(tconst);
      satvp = esat(p72.ptemp(:,indTime));

      %% see HumidityMeasures.pdf     e = q p/(0.622 + 0.378 q)    
      vp = p72.gas_1(:,indTime) .* p72.plevs(:,indTime) ./(0.622 + 0.378* p72.gas_1(:,indTime));         

      figure(1); pcolor(rh);             colorbar; colormap(jet);                      shading interp; title('RH direct from p,T,q');
      figure(2); pcolor(vp./satvp);      colorbar; colormap(jet);                      shading interp; title('RH from vp and satvp');
      figure(3); pcolor(rh - vp./satvp); colorbar; colormap(usa2); caxis([-1 +1]*0.5); shading interp; title('RH from (vp and satvp) compared to direct');

      % so satvp = qsat p/(0.622 + 0.378 qsat)
      % so qsat = 0.622 satvp / (p - 0.378 satvp)
      qsat = 0.622 * satvp ./ (p72.plevs(:,indTime) - 0.378 * satvp);

      p72.gas_1(:,indTime) = wvconst + dq + ((iTime-1)*RHincrease/12)*qsat;
    end   

    iPlot = -1;
    if iPlot > 0
      figure(1); pcolor(tconst);                       shading interp; colorbar; title('has T changed?');   colormap(jet); set(gca,'ydir','reverse')
      figure(1); pcolor(tconst-p72.ptemp(:,indTime));  shading interp; colorbar; title('has T changed?');   colormap(usa2); caxis([-1 +1]*10); set(gca,'ydir','reverse')
      figure(2); pcolor(wvconst);                      shading interp; colorbar; title('has WV changed?');   colormap(jet); set(gca,'ydir','reverse')
      figure(2); pcolor(wvconst-p72.gas_1(:,indTime)); shading interp; colorbar; title('has WV changed?');   colormap(usa2); caxis([-1 +1]*1e-3); set(gca,'ydir','reverse')
      figure(3); pcolor(wvconst-initialWV(:,indTime)); shading interp; colorbar; title('has WV changed?');   colormap(usa2); caxis([-1 +1]*1e-3); set(gca,'ydir','reverse')
      figure(4); pcolor(rhconst-rhnew); colorbar; title('has RH changed?');  colormap(usa2); caxis([-1 +1]/10); shading interp; set(gca,'ydir','reverse')
      pause(0.1)
    end
    rhnew = convert_humidity (p72.plevs(:,indTime),p72.ptemp(:,indTime),p72.gas_1(:,indTime),'specific humidity','relative humidity');
    drhsave(iTime) = rhnew(20,36)-rhconst(20,36);
    figure(6); clf; plot(drhsave); title(['idRH = ' num2str(idRH) ' \delta RH at center of 37 x 72 set']); set(gca,'fontsize',10)
    figure(5); clf; pcolor(rhnew-rhconst); caxis([0 0.25]); shading interp; set(gca,'ydir','reverse'); colormap(usa2); colorbar; 
      title(['idRH = ' num2str(idRH) ' ... RHnew-RH0 at time ' num2str(iTime)])
      %disp('ret to continue'); pause

      pause(0.1)
  end
  p72.gas_1 = max(0,p72.gas_1);
elseif idRH == 5
  disp('idRH == 5, use clouds and everything')
  p72 = fix_clouds_as_needed(p72);
elseif idRH == 7
  disp('idRH == 7, do nothing and had no clouds')
end
