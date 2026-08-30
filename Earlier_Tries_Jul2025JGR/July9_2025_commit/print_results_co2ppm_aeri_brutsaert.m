data = [p5.rlon; p5.rlat; p5.salti; p5.landfrac; p5.spres; p5.stemp; p5.mmw; p5.olr_clr; -p5.ilr_clr+5.67e-8*p5.stemp.^4; p5.d2m;  p5.t2m; ppmvAVG; ecRad0.clr; ecRadWV; ecRadCO2; ecRadT; ecRadtotal];
[yy,mm,dd] = tai2utcSergio(nanmean(p5.rtime));

if rCO2PPM == -1
  fid = fopen(['akleidon_' num2str(iMonthSoFar,'%03d') '.txt'],'w');
else
  fid = fopen(['akleidon_' num2str(iMonthSoFar,'%03d') '_co2ppm_' num2str(rCO2PPM) '.txt'],'w');
end

fprintf(fid,'saving data for yy/mm/dd = %4i/%2i/%2i \n',yy,mm,dd);
fprintf(fid,'   rlon      rlat     surf_alt   landfrac   surf_pres   surf_temp  mm_water   TOA_OLR_ERA5    GND_ILR_ERA5  Tdew2m   Tair2m   | CO2amt    | Rld_ecRad   Rld_WV_ecRad  Rld_CO2_ecRad  Rld_Tz_ecRad   Rld_Total_ecRad \n');
fprintf(fid,'  [deg]      [deg]      [m]        [0-1]      [mb]        [K]        [mm]        [W/m2]          [W/m2]      [K]       [K]    |  [ppmv]   |   [W/m2]       [W/m2]        [W/m2]         [W/m2]           [W/m2]  \n');
fprintf(fid,'    1          2         3          4           5          6          7            8               9         10        11         12           13            14           15              16               17    \n');
fprintf(fid,'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%8.3f  %8.3f    %8.3f    %5.4f    %8.3f    %8.3f   %8.3f      %8.3f    %8.3f    %8.3f  %8.3f   |  %8.3f |  %8.3f   %8.3f      %8.3f       %8.3f        %8.3f  \n',data); 
fclose(fid);

figure(1); clf; scatter_coast(data(1,:),data(2,:),100,data(4,:)); colormap jet; title('Landfrac')
figure(2); clf; plot(data(2,:),data(14,:),'b',data(2,:),data(15,:),'g',data(2,:),data(16,:),'r',data(2,:),data(17,:),'k','linewidth',2);
  plotaxis2; xlabel('Latitude'); ylabel('Flux Change W/m2');
  hl = legend('WV','CO2','T','Total','location','best','fontsize',10);
xlim([-max(abs(p5.rlat)) +max(abs(p5.rlat))])

%{
 a300 = load('akleidon_120_co2ppm_300.txt');
 a400 = load('akleidon_120.txt');
 plot(a300(:,12),a300(:,13)-a400(:,13))
 plot(a300(:,1),a300(:,13)-a400(:,13))
 plot(a300(:,2),a300(:,13)-a400(:,13))
 plot(a300(:,2),a300(:,13)-a400(:,13)); xlabel('Latitude'); ylabel('Rld(300 ppm) - Rld(400 ppm)  W/m2')
 plot(a300(:,2),a300(:,13)-a400(:,13)); xlabel('Latitude'); title('Rld(300 ppm) - Rld(400 ppm)'); ylabel('\delta (Rld) W/m2')
 plotaxis2;
%}

