function y = plot_tracegas_VS_latbin_anomaly_runs(file1,file2);

%{
f9a = 'anomaly_0dayavg_results_strow_latbin_unc_V3co2n2och4jac_afglprofs.mat';     %% strow jacs, had a big bug in jacs
f9b = 'anomaly_0dayavg_cal_results_strow_latbin_unc_V3co2n2och4jac_afglprofs.mat'; %% strow jacs, had a big bug in jacs

f10a = 'anomaly_0dayavg_results_strow_latbin_unc_afglprofs.mat';     %% time varying unc, sergio jacs
f10b = 'anomaly_0dayavg_cal_results_strow_latbin_unc_afglprofs.mat'; %% time varying unc, sergio jacs

f11a = 'anomaly_0dayavg_results_strow_latbin_unc_afglprofs_V3co2n2och4jac_fixedjac.mat';     %% strow jacs, fixed jac bugs
f11b = 'anomaly_0dayavg_cal_results_strow_latbin_unc_V3co2n2och4jac_afglprofs_fixedjac.mat'; %% strow jacs, fixed jac bugs

x9  = plot_tracegas_VS_latbin_anomaly_runs(f9a,f9b);
x10 = plot_tracegas_VS_latbin_anomaly_runs(f10a,f10b);
x11 = plot_tracegas_VS_latbin_anomaly_runs(f11a,f11b);

figure(1); plot(x9.iaLat,x9.stemp,'bo-',x10.iaLat,x10.stemp,'rx-','linewidth',2); grid; axis([-60 +60 -0.03 +0.03])
figure(2); plot(x9.iaLat,x9.co2,'bo-',x10.iaLat,x10.co2,'rx-','linewidth',2); grid; axis([-60 +60 1.75 2.45])
figure(3); plot(x9.iaLat,x9.n2o,'bo-',x10.iaLat,x10.n2o,'rx-','linewidth',2); grid; axis([-60 +60 0.75 1.5])
figure(4); plot(x9.iaLat,x9.ch4,'bo-',x10.iaLat,x10.ch4,'rx-','linewidth',2); grid; axis([-60 +60 4.00 8.00])
figure(5); plot(x9.iaLat,x9.co2save_lat_t(:,1:3)-x10.co2save_lat_t(:,1:3),'bx-',...
                x9.iaLat,x9.co2save_lat_t(:,4:6)-x10.co2save_lat_t(:,4:6),'ro-','linewidth',2); grid; axis([-60 +60 -5 +5])
  title('(A-B) (b) first3 times (r) last 3 times')
figure(6); plot(x9.iaLat,x9.co2save_lat_t(:,4:6)-x9.co2save_lat_t(:,1:3),'bx-',...
                x9.iaLat,x10.co2save_lat_t(:,4:6)-x10.co2save_lat_t(:,1:3),'ro-','linewidth',2); grid; axis([-60 +60 0 40])
  title('( (b) A last-first 3 times (r) B last-first 3 times')

figure(1); plot(x11.iaLat,x11.stemp,'bo-',x10.iaLat,x10.stemp,'rx-','linewidth',2); grid; axis([-60 +60 -0.03 +0.03])
figure(2); plot(x11.iaLat,x11.co2,'bo-',x10.iaLat,x10.co2,'rx-','linewidth',2); grid; axis([-60 +60 1.75 2.45])
figure(3); plot(x11.iaLat,x11.n2o,'bo-',x10.iaLat,x10.n2o,'rx-','linewidth',2); grid; axis([-60 +60 0.75 1.5])
figure(4); plot(x11.iaLat,x11.ch4,'bo-',x10.iaLat,x10.ch4,'rx-','linewidth',2); grid; axis([-60 +60 4.00 8.00])
figure(5); plot(x11.iaLat,x11.co2save_lat_t(:,1:3)-x10.co2save_lat_t(:,1:3),'bx-',...
                x11.iaLat,x11.co2save_lat_t(:,4:6)-x10.co2save_lat_t(:,4:6),'ro-','linewidth',2); grid; axis([-60 +60 -5 +5])
  title('(A-B) (b) first3 times (r) last 3 times')
figure(6); plot(x11.iaLat,x11.co2save_lat_t(:,4:6)-x11.co2save_lat_t(:,1:3),'bx-',...
                x11.iaLat,x10.co2save_lat_t(:,4:6)-x10.co2save_lat_t(:,1:3),'ro-','linewidth',2); grid; axis([-60 +60 0 40])
  title('( (b) A last-first 3 times (r) B last-first 3 times')
%}

%% look at compare_anomaly_runs
%% file1 usually OBS and file2 CAL so file3 = difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSmooth = input('Enter how many timesteps to smooth over (each timestep = 16 days = 0.0438 yr; default = 11 timesteps = 0.5 yr) : ');

a1 = load(file1);
a2 = load(file2); 

addpath /home/sergio/MATLABCODE/PLOTTER
iaLat = equal_area_spherical_bands(20);
iaLat = 0.5 * (iaLat(1:end-1) + iaLat(2:end));
tropics = find(abs(iaLat) < 30);

if nargin == 2
  %iYes = input('pretend a3 = a1-a2 (-1/+1) ? ');
  iYes = 1;
  if iYes > 0
    a3.stemp = a1.stemp - a2.stemp;
    a3.co2 = a1.co2 - a2.co2;
    a3.n2o = a1.n2o - a2.n2o;
    a3.ch4 = a1.ch4 - a2.ch4;

    if ~isfield(a1,'topts')
      a1.topts.nloop = -9999;
    end
    if ~isfield(a2,'topts')
      a2.topts.nloop = -9999;
    end

    if ~isfield(a1,'topts')
      a1.topts.nloop = -9999;
    end
    if ~isfield(a2,'topts')
      a2.topts.nloop = -9999;
    end

    if isfield(a1,'topts') & isfield(a2,'topts')
      if isfield(a1.topts,'nloop') & isfield(a2.topts,'nloop')
        a3.topts.nloop = min(a1.topts.nloop,a2.topts.nloop);
      else
        a3.topts.nloop = -9999;
      end
    else
      a3.topts.nloop = -9999;
    end
    a3.okdates = a2.okdates;
  end
end


iOK = +1;
if iOK > 0
  junk = [sum(sum(a1.co2(tropics,:)-a2.co2(tropics,:))) sum(sum(a1.stemp(tropics,:)-a2.stemp(tropics,:)))];
  fprintf(1,'sum(diff 1,2 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))
  junk = [sum(sum(a1.co2(tropics,:)-a3.co2(tropics,:))) sum(sum(a1.stemp(tropics,:)-a3.stemp(tropics,:)))];
  fprintf(1,'sum(diff 1,3 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))
  junk = [sum(sum(a2.co2(tropics,:)-a3.co2(tropics,:))) sum(sum(a2.stemp(tropics,:)-a3.stemp(tropics,:)))];
  fprintf(1,'sum(diff 2,3 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))

  fprintf(1,'co2   at timestep 364, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.co2(21,364) a2.co2(21,364) a3.co2(21,364)])
  fprintf(1,'stemp at timestep 364, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.stemp(21,364) a2.stemp(21,364) a3.stemp(21,364)])
  fprintf(1,'nloop at timestep 364, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.topts.nloop a2.topts.nloop a3.topts.nloop])

  figure(1); plot(a1.okdates,smooth(nanmean(a1.stemp(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.stemp(tropics,:)),iSmooth),'g.-',...
                  a1.okdates,smooth(nanmean(a3.stemp(tropics,:)),iSmooth),'r.-'); title('tropical stemp');
  figure(2); plot(a1.okdates,smooth(nanmean(a1.co2(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.co2(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.co2(tropics,:)),iSmooth),'r.-'); title('tropical co2');
  figure(3); plot(a1.okdates,smooth(nanmean(a1.ch4(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.ch4(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.ch4(tropics,:)),iSmooth),'r.-'); title('tropical ch4');
  figure(4); plot(a1.okdates,smooth(nanmean(a1.n2o(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.n2o(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.n2o(tropics,:)),iSmooth),'r.-'); title('tropical n2o');

  disp(' ')
  disp('Tropical avg rates')

  P1 = polyfit(a2.okdates',smooth(nanmean(a1.stemp(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.stemp(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.stemp(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed stemp rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.co2(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.co2(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.co2(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed co2 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.n2o(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.n2o(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.n2o(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed n2o rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.ch4(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.ch4(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.ch4(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed ch4 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.stemp(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.stemp(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.stemp(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed stemp rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.co2(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.co2(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.co2(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed co2 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.n2o(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.n2o(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.n2o(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed n2o rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.ch4(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.ch4(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.ch4(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed ch4 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  for ii = 1 : 4
    figure(ii); ax = axis; ax(1) = min(a1.okdates); ax(2) = max(a1.okdates); axis(ax); grid on;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ll = 1 : 40
    ind = 3 : length(a1.okdates)-2;

    gah = smooth(a1.stemp(ll,:),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a2.stemp(ll,:),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a3.stemp(ll,:),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
    stemp_1(ll) = P1(1);
    stemp_2(ll) = P2(1);
    stemp_3(ll) = P3(1);

    gah = smooth(a1.co2(ll,:),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a2.co2(ll,:),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a3.co2(ll,:),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
    co2_1(ll) = P1(1);
    co2_2(ll) = P2(1);
    co2_3(ll) = P3(1);

    gah = smooth(a1.n2o(ll,:),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a2.n2o(ll,:),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a3.n2o(ll,:),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
    n2o_1(ll) = P1(1);
    n2o_2(ll) = P2(1);
    n2o_3(ll) = P3(1);

    gah = smooth(a1.ch4(ll,:),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a2.ch4(ll,:),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
    gah = smooth(a3.ch4(ll,:),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
    ch4_1(ll) = P1(1);
    ch4_2(ll) = P2(1);
    ch4_3(ll) = P3(1);
  end

  figure(5); plot(iaLat,stemp_3,'x-','linewidth',2); title('stemp rates K/yr'); axis([-60 +60 -0.02 +0.02]); grid
  figure(6); plot(iaLat,co2_3,'x-','linewidth',2); title('co2 rates ppm/yr'); axis([-60 +60 1.8 2.2]); grid
  figure(7); plot(iaLat,n2o_3,'x-','linewidth',2); title('n2o rates ppb/yr'); axis([-60 +60 0.5 1.5]); grid
  figure(8); plot(iaLat,ch4_3,'x-','linewidth',2); title('ch4 rates ppb/yr'); axis([-60 +60 5 7]); grid

  figure(5); plot(iaLat,smooth(stemp_3,5),'x-','linewidth',2); title('stemp rates K/yr'); axis([-60 +60 -0.02 +0.02]); grid
  figure(6); plot(iaLat,smooth(co2_3,5),'x-','linewidth',2); title('co2 rates ppm/yr'); axis([-60 +60 1.8 2.2]); grid
  figure(7); plot(iaLat,smooth(n2o_3,5),'x-','linewidth',2); title('n2o rates ppb/yr'); axis([-60 +60 0.5 1.5]); grid
  figure(8); plot(iaLat,smooth(ch4_3,5),'x-','linewidth',2); title('ch4 rates ppb/yr'); axis([-60 +60 5 7]); grid

  y.iaLat = iaLat;
  y.stemp = smooth(stemp_3,5);
  y.co2 = smooth(co2_3,5);
  y.n2o = smooth(n2o_3,5);
  y.ch4 = smooth(ch4_3,5);
  y.co2save_lat_t = a3.co2(:,[3 4 5 361 362 363]);

  figure(1); plot(iaLat,a3.co2(:,[3 4 5 361 362 363])); grid; axis([-60 +60 -10 +60]); 
end

disp(' ')
