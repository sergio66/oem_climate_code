addpath /home/sergio/MATLABCODE/PLOTTER
iOffSet = 10;
offT = find(A.okdates >= C.okdates(1),1);

figure(iOffSet+1); clf
offY  = A.co2(offT);
offY1 = A.co2_1(offT);
  h1 = subplot(211); 
    plot(C.okdates,C.co2+offY,'b',A.okdates,A.co2,'r','linewidth',1);
    hl = legend('CRIS obs-cal','AIRS obs-cal','location','best'); set(hl,'fontsize',8);
    ylabel('CO2'); grid
  h2 = subplot(212); 
    plot(C.okdates,C.co2_1+offY1,'b',A.okdates,A.co2_1,'r','linewidth',1);
    hl = legend('CRIS obs','AIRS obs','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('CO2'); grid
  adjust21(h1,h2,'even')

figure(iOffSet+2); clf
offY  = A.ch4(offT);
offY1 = A.ch4_1(offT);
  h1 = subplot(211); 
    plot(C.okdates,C.ch4+offY,'b',A.okdates,A.ch4,'r','linewidth',1);
    hl = legend('CRIS obs-cal','AIRS obs-cal','location','best'); set(hl,'fontsize',8);
    ylabel('CH4'); grid
  h2 = subplot(212); 
    plot(C.okdates,C.ch4_1+offY1,'b',A.okdates,A.ch4_1,'r','linewidth',1);
    hl = legend('CRIS obs','AIRS obs','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('CH4'); grid
  adjust21(h1,h2,'even')

figure(iOffSet+3); clf
offY  = A.stemp(offT);
offY1 = A.stemp_1(offT);
offY  = 0;
offY1 = 0;
  h1 = subplot(211); 
    plot(C.okdates,C.stemp+offY,'b',A.okdates,A.stemp,'r','linewidth',1);
    hl = legend('CRIS obs-cal','AIRS obs-cal','location','best'); set(hl,'fontsize',8);
    ylabel('STEMP'); grid
  h2 = subplot(212); 
    plot(C.okdates,C.stemp_1+offY1,'b',A.okdates,A.stemp_1,'r','linewidth',1);
    hl = legend('CRIS obs','AIRS obs','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('STEMP'); grid
  adjust21(h1,h2,'even')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ../AIRS_new_clear_scan_August2019_AMT2020PAPER/f2645.mat
load ../CRIS_new_clear_scan_January2020/f1305.mat

iA791 = find(f2645 >= 791,1);
iA792 = find(f2645 >= 792,1);
iC791 = find(f1305 >= 791,1);
iC792 = find(f1305 >= 792,1);
iC791 = find(f1305 >= 791-0.5,1);
iC792 = find(f1305 >= 792-0.5,1);

offY  = btA.raaObs(iA791,offT);
offY1 = btA.raaCal(iA791,offT);
offY  = 0;
offY1 = 0;
figure(iOffSet+4); clf
  h1 = subplot(211); 
    plot(C.okdates,btC.raaObs(iC791,:)+offY,'b',A.okdates,btA.raaObs(iA791,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best'); set(hl,'fontsize',8);
    ylabel('BT791 obs'); grid
  h2 = subplot(212); 
    plot(C.okdates,btC.raaCal(iC791,:)+offY1,'b',A.okdates,btA.raaCal(iA791,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('BT791 cal'); grid
  adjust21(h1,h2,'even')

figure(iOffSet+5); clf
  h1 = subplot(211); 
    plot(C.okdates,btC.raaObs(iC791-1,:)-btC.raaObs(iC792+1,:)+0.4,'b',A.okdates,btA.raaObs(iA791,:)-btA.raaObs(iA792,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best'); set(hl,'fontsize',8);
    ylabel('BT791-792 obs'); grid
  h2 = subplot(212); 
    plot(C.okdates,btC.raaCal(iC791-1,:)-btC.raaCal(iC792+1,:)+0.4,'b',A.okdates,btA.raaCal(iA791,:)-btA.raaCal(iA792,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('BT791-792 cal'); grid
  adjust21(h1,h2,'even')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(iOffSet+6); clf
iA1305 = find(f2645 >= 1305,1);
iC1305 = find(f1305 >= 1305,1);

offY  = btA.raaObs(iA1305,offT);
offY1 = btA.raaCal(iA1305,offT);
offY  = 0;
offY1 = 0;
  h1 = subplot(211); 
    plot(C.okdates,btC.raaObs(iC1305,:)+offY,'b',A.okdates,btA.raaObs(iA1305,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best'); set(hl,'fontsize',8);
    ylabel('BT1305 obs'); grid
  h2 = subplot(212); 
    plot(C.okdates,btC.raaCal(iC1305,:)+offY1,'b',A.okdates,btA.raaCal(iA1305,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('BT1305 cal'); grid
  adjust21(h1,h2,'even')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(iOffSet+7); clf

iA1231 = find(f2645 >= 1231,1);
iC1231 = find(f1305 >= 1231,1);

offY  = btA.raaObs(iA1231,offT);
offY1 = btA.raaCal(iA1231,offT);
offY  = 0;
offY1 = 0;
  h1 = subplot(211); 
    plot(C.okdates,btC.raaObs(iC1231,:)+offY,'b',A.okdates,btA.raaObs(iA1231,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best'); set(hl,'fontsize',8);
    ylabel('BT1231 obs'); grid
  h2 = subplot(212); 
    plot(C.okdates,btC.raaCal(iC1231,:)+offY1,'b',A.okdates,btA.raaCal(iA1231,:),'r','linewidth',1);
    hl = legend('CRIS','AIRS','location','best');  set(hl,'fontsize',8);
    xlabel('time'); ylabel('BT1231 cal'); grid
  adjust21(h1,h2,'even')


