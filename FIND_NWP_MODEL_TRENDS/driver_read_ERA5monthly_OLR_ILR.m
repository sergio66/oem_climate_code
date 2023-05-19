addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/

noPi = 3.1415; %% should we have this factor in the solar bean???
noPi = 1.0;   %% I think it is solar beam, so no solid angle integral

fprintf(1,'stefan boltzmann law at 255 K = %8.6f W/m2 \n',5.6704e-8 * 255^4)

%% 6.785087652174316e-5 = pi(sun diam/sun dist)/asl/rta/kcarta_sergio2 = pi*((0.6951e9/149.57e9)/asl/rta/kcarta_sergio2)
rOmegaSun = 6.785087652174316e-5; %% version on /home/sergio/KCARTA/INCLUDE/kcarta.param
rSWalbedo = 0.3;                  %% refelectance from clouds

x = 1:0.1:3000; figure(49); y = ttorad(x,285); plot(x,y); xlabel('Wavenumber cm-1'); ylabel('Rad mW/m2/sr/cm-1')
x = 1:0.1:3000; figure(49); y = ttorad(x,255); plot(x,y); xlabel('Wavenumber cm-1'); ylabel('Rad mW/m2/sr/cm-1')
  fprintf(1,'integrating upwelling thermal planck gives %8.6f W/m2 \n',sum(y)*mean(diff(x))/1000*pi)

%xs = [0.1 0.9]; xs = sort(10000./xs); xs = linspace(xs(1),xs(2),10000); figure(50); 
%  ys = ttorad(xs,5800); plot(10000./xs,ys*rOmegaSun*(1-rSWalbedo)); sum(ys)*abs(mean(diff(xs)))*rOmegaSun/1000*noPi
xs = [0.1 10];   xs = sort(10000./xs); xs = linspace(xs(1),xs(2),10000); figure(50); 
  ys = ttorad(xs,5800); plot(10000./xs,ys*rOmegaSun*(1-rSWalbedo)); 
  fprintf(1,'integrating incoming (wavenumber,Plank(wavenumber,5800))*rOmegaSun gives %8.6f W/m2 \n',sum(ys)*abs(mean(diff(xs)))*rOmegaSun/1000*noPi)

% w      is in cm-1, rin  is in mW/m2/cm-1/sr
% lambda is in um,   rout is in W /m2/um/sr
[lambda,rout] = rads_wnum2lamda(xs,ys); figure(51); 
  plot(lambda,rout*rOmegaSun*(1-rSWalbedo)*noPi);
  plot(lambda,rout*rOmegaSun*noPi); xlim([0 4]); %% at TOA; 
  xlabel('Wavelength um'); ylabel('Rad W/m2/sr/um')
  [Y,I] = sort(lambda);
  fprintf(1,'integrating incoming (lambda,Plank(lambda,5800))*rOmegaSun with trapz  gives %8.6f W/m2 \n',trapz(lambda(I),rout(I))*rOmegaSun*noPi)

figure(52); plot(x,y,xs,ys*rOmegaSun); xlabel('Wavenumber cm-1'); ylabel('Rad mW/m2/sr/cm-1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%x = read_netcdf_lls('/asl/models/era5_avg/INCOMING/2022-03_lev.nc');  
x2 = read_netcdf_lls('/asl/models/era5_avg/INCOMING/2022-03_sfc.nc');

%x = read_netcdf_lls('/asl/models/era5_avg/INCOMING/2022-10_lev.nc');  
x2 = read_netcdf_lls('/asl/models/era5_avg/INCOMING/2022-10_sfc.nc');

[Y,X] = meshgrid(x2.latitude,x2.longitude);

jett = jet; jett(1,:) = 1;

figure(40); colormap(jet);  Z = squeeze(x2.mtnlwrfcs(:,:,6)); simplemap(Y,X,abs(Z)); title('LW OLR clear sky')
figure(41); colormap(jet);  Z = squeeze(x2.mtnlwrf(:,:,6));   simplemap(Y,X,abs(Z)); title('LW OLR allsky')
figure(42); colormap(jet);  Z = squeeze(x2.msnlwrfcs(:,:,6)); simplemap(Y,X,abs(Z)); title('LW ILR clear sky')
figure(43); colormap(jet);  Z = squeeze(x2.msnlwrf(:,:,6));   simplemap(Y,X,abs(Z)); title('LW ILR allsky')
figure(44); colormap(jet);  Z = squeeze(x2.skt(:,:,5));       simplemap(Y,X,Z); title('SKT (K)')
figure(45); colormap(jett); Z = squeeze(x2.mtnswrfcs(:,:,6)); simplemap(Y,X,abs(Z)); title('SW up TOA clear sky')
figure(46); colormap(jett); Z = squeeze(x2.mtnswrf(:,:,6));   simplemap(Y,X,abs(Z)); title('SW up TOA allsky')
figure(47); colormap(jett); Z = squeeze(x2.msnswrfcs(:,:,6)); simplemap(Y,X,abs(Z)); title('SW dn GND clear sky')
figure(48); colormap(jett); Z = squeeze(x2.msnswrf(:,:,6));   simplemap(Y,X,abs(Z)); title('SW dn GND allsky')

%{
https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings
mean_surface_downward_long_wave_radiation_flux		     msdwlwrf
mean_surface_net_long_wave_radiation_flux		     msnlwrf     Y
mean_top_net_long_wave_radiation_flux			     mtnlwrf     Y
mean_top_net_long_wave_radiation_flux_clear_sky		     mtnlwrfcs   Y
mean_surface_net_long_wave_radiation_flux_clear_sky	     msnlwrfcs   Y
mean_surface_downward_long_wave_radiation_flux_clear_sky     msdwlwrfcs
%}
