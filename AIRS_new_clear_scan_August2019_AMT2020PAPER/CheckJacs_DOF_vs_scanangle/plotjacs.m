
figure(4); 
  plot(fKc(1:2378),sum(WV_45ang(1:2378,:)'),'b',fKc(1:2378),sum(WV_22ang(1:2378,:)'),'g',fKc(1:2378),sum(WV_00ang(1:2378,:)'),'r',...
       fKc(1:2378),sum(WV_avg(1:2378,:)'),'kx-');
  hl = legend('45 deg','22 deg','0 deg','avg','location','best'); title('col WV')
figure(5); 
  plot(fKc(1:2378),sum(O3_45ang(1:2378,:)'),'b',fKc(1:2378),sum(O3_22ang(1:2378,:)'),'g',fKc(1:2378),sum(O3_00ang(1:2378,:)'),'r',...
       fKc(1:2378),sum(O3_avg(1:2378,:)'),'kx-');
  hl = legend('45 deg','22 deg','0 deg','avg','location','best'); title('col O3')
figure(6); 
  plot(fKc(1:2378),sum(T_45ang(1:2378,:)'),'b',fKc(1:2378),sum(T_22ang(1:2378,:)'),'g',fKc(1:2378),sum(T_00ang(1:2378,:)'),'r',...
       fKc(1:2378),sum(T_avg(1:2378,:)'),'kx-');
  hl = legend('45 deg','22 deg','0 deg','avg','location','best'); title('col T')
figure(7); 
  plot(fKc(1:2378),stemp_45ang(1:2378),'b',fKc(1:2378),stemp_22ang(1:2378),'g',fKc(1:2378),stemp_00ang(1:2378),'r',...
       fKc(1:2378),stemp_avg(1:2378),'kx-');
  hl = legend('45 deg','22 deg','0 deg','avg','location','best'); title('Surf Temp')

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); 
  plot(fKc(1:2378),sum(WV_22ang(1:2378,:)'),'b',fKc(1:2378),sum(WV_avg(1:2378,:)'),'r');
  hl = legend('22 deg','avg','location','best'); title('col WV')
figure(5); 
  plot(fKc(1:2378),sum(O3_22ang(1:2378,:)'),'b',fKc(1:2378),sum(O3_avg(1:2378,:)'),'r');
  hl = legend('22 deg','avg','location','best'); title('col O3')
figure(6); 
  plot(fKc(1:2378),sum(T_22ang(1:2378,:)'),'b',fKc(1:2378),sum(T_avg(1:2378,:)'),'r');
  hl = legend('22 deg','avg','location','best'); title('col T')
figure(7); 
  plot(fKc(1:2378),stemp_22ang(1:2378),'b',fKc(1:2378),stemp_avg(1:2378),'r');
  hl = legend('22 deg','avg','location','best'); title('Sutf Temp')

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); 
  plot(fKc(1:2378),sum(WV_45ang(1:2378,:)'),'b',fKc(1:2378),sum(WV_22ang(1:2378,:)'),'g',fKc(1:2378),sum(WV_00ang(1:2378,:)'),'r')
  hl = legend('45 deg','22 deg','0 deg','location','best'); title('kCARTA col WV')
figure(5); 
  plot(fKc(1:2378),sum(O3_45ang(1:2378,:)'),'b',fKc(1:2378),sum(O3_22ang(1:2378,:)'),'g',fKc(1:2378),sum(O3_00ang(1:2378,:)'),'r')
  hl = legend('45 deg','22 deg','0 deg','location','best'); title('kCARTA col O3')
figure(6); 
  plot(fKc(1:2378),sum(T_45ang(1:2378,:)'),'b',fKc(1:2378),sum(T_22ang(1:2378,:)'),'g',fKc(1:2378),sum(T_00ang(1:2378,:)'),'r')
  hl = legend('45 deg','22 deg','0 deg','location','best'); title('kCARTA col T')
figure(7); 
  plot(fKc(1:2378),stemp_45ang(1:2378),'b',fKc(1:2378),stemp_22ang(1:2378),'g',fKc(1:2378),stemp_00ang(1:2378),'r')
  hl = legend('45 deg','22 deg','0 deg','location','best'); title('kCARTA Surf Temp')

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
