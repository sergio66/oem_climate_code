%{
https://en.wikipedia.org/wiki/Broyden%27s_method

J(n) = J(n-1) +  dy(n) - J(n-1) dx(n)
                 -------------------- dx(n)'
                       |dx(n)|^2

J(n) = J(n-1) +  dy(n) - J(n-1) dx(n)
                 -------------------- dx(n)'
                     |dx(n)'*dx(n)|

units dy = BT
      J = d(BT)/dx
thus  BT - Bt/dx dx  dx      = dBT/dx = jacobian
      -------------
          dx^2

also look at /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/jim_lambers_numericalanalysis.pdf pg 323
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/h4tools

%[h,ha,pbefore,pa] = rtpread('pbeforeavg.op.rtp');  %% before mar 30, 2020
%[h,ha,pafter,pa]  = rtpread('pafteravg.op.rtp');   %% before mar 30, 2020

[h,ha,pbefore,pa] = rtpread('pbeforeavg2.op.rtp');  %% after apr 23, 2020
[h,ha,pafter,pa]  = rtpread('pafteravg2.op.rtp');   %% after apr 23, 2020

iiBin = 30;
iiBin = 20;

w0 = load(['PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);    w1 = load(['PAFTER/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);
a0 = load(['PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']); a1 = load(['PAFTER/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']);

figure(1); plot(w0.f,sum(w0.jacWV'),'b',w0.f,sum(w1.jacWV'-w0.jacWV'),'r'); title('WV jac')
figure(2); plot(w0.f,sum(w0.jacT'),'b',w0.f,sum(w1.jacT'-w0.jacT'),'r');    title('T jac');
figure(3); plot(a0.f,sum(a0.jacO3'),'b',a0.f,sum(a1.jacO3'-a0.jacO3'),'r'); title('O3 jac');
figure(4); plot(a0.f,a0.jacST,'b',a0.f,a1.jacST-a0.jacST,'r');              title('STEMP jac')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('look at save_the_rtp, "pafteravg", for zonal plots about the mean profile differences ....')
disp('look at save_the_rtp, "pafteravg", for zonal plots about the mean profile differences ....')
disp('look at save_the_rtp, "pafteravg", for zonal plots about the mean profile differences ....')

%{
f = w0.f;
colWVjac0  = sum(w0.jacWV');  colWVjac1  = sum(w1.jacWV'); 
colTjac0   = sum(w0.jacT');   colTjac1   = sum(w1.jacT'); 
colO3jac0  = sum(a0.jacO3');  colO3jac1  = sum(a1.jacO3'); 
colCO2jac0 = (a0.jacCO2'); colCO2jac1 = (a1.jacCO2'); 
colCH4jac0 = (a0.jacCH4'); colCH4jac1 = (a1.jacCH4'); 
colN2Ojac0 = (a0.jacN2O'); colN2Ojac1 = (a1.jacN2O'); 
stjac0     = a0.jacST;        stjac1    = a1.jacST;
[h,psave0] = subset_rtp(h,pbefore,[],[],iiBin);
[h,psave1] = subset_rtp(h,pafter,[],[],iiBin);
  psave0.std_stemp = pbefore.std_stemp(iiBin);
  psave0.std_ptemp = pbefore.std_ptemp(:,iiBin);
  psave0.std_gas_1 = pbefore.std_gas_1(:,iiBin);
  psave0.std_gas_3 = pbefore.std_gas_3(:,iiBin);
  psave0.std_stemp = pafter.std_stemp(iiBin);
  psave0.std_ptemp = pafter.std_ptemp(:,iiBin);
  psave0.std_gas_1 = pafter.std_gas_1(:,iiBin);
  psave0.std_gas_3 = pafter.std_gas_3(:,iiBin);
%save for_amt2020_review.mat f col* psave0 psave1 iiBin     %% before Mar 30, 2020
save for_amt2020_review2.mat f col* psave0 psave1 iiBin     %% after  Apr 23, 2020
%}
error(';ksjgksjg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(w0.f >= 800 & w0.f <= 980);
bt0 = rad2bt(w0.f,w0.rad);  bt1 = rad2bt(w0.f,w1.rad);  dy = bt1-bt0; 
dx = pafter.stemp(iiBin)-pbefore.stemp(iiBin);
dJall = (dy - a0.jacST * dx)/(dx*dx) * dx;
figure(5); plot(w0.f,a1.jacST-a0.jacST,'r',w0.f,dJall,'k'); title('Stemp')
%dJ = (dy(ind) - a0.jacST(ind) * dx)/(abs(dx)^2) * dx;
%figure(5); plot(w0.f(ind),a1.jacST(ind)-a0.jacST(ind),'r',w0.f(ind),dJ,'kx-',w0.f,dJall,'g'); 

bt0 = rad2bt(w0.f,w0.rad);  bt1 = rad2bt(w0.f,w1.rad);  dy = bt1-bt0; 
dx = pafter.ptemp(1:97,iiBin)-pbefore.ptemp(1:97,iiBin); dx = flipud(dx);
dJ = (dy - a0.jacT * dx)/(dx'*dx) * dx';
figure(6); plot(w0.f,sum(a1.jacT'-a0.jacT'),'r',w0.f,sum(dJ'),'k'); title('T')

bt0 = rad2bt(w0.f,w0.rad);  bt1 = rad2bt(w0.f,w1.rad);  dy = bt1-bt0; 
dx = log10(pafter.gas_1(1:97,iiBin))-log10(pbefore.ptemp(1:97,iiBin)); dx = flipud(dx);
dJ = (dy - w0.jacWV * dx)/(dx'*dx) * dx';
figure(7); plot(w0.f,sum(w1.jacWV'-w0.jacWV'),'r',w0.f,sum(dJ'),'k'); title('WV')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% da big one

bt0 = rad2bt(w0.f,w0.rad);  bt1 = rad2bt(w0.f,w1.rad);  dy = bt1-bt0; 
X0 = [pbefore.stemp(iiBin); pbefore.ptemp(1:97,iiBin); log(pbefore.gas_1(1:97,iiBin)); log(pbefore.gas_3(1:97,iiBin))];
X1 = [pafter.stemp(iiBin);  pafter.ptemp(1:97,iiBin);  log(pafter.gas_1(1:97,iiBin));  log(pafter.gas_3(1:97,iiBin))];
X0 = [pbefore.stemp(iiBin); flipud(pbefore.ptemp(1:97,iiBin)); log(flipud(pbefore.gas_1(1:97,iiBin))); log(flipud(pbefore.gas_3(1:97,iiBin))) ];
X1 = [pafter.stemp(iiBin);  flipud(pafter.ptemp(1:97,iiBin));  log(flipud(pafter.gas_1(1:97,iiBin)));  log(flipud(pafter.gas_3(1:97,iiBin)))] ;
dx = X1-X0;
jac0 = [a0.jacST a0.jacT w0.jacWV a0.jacO3];
jac1 = [a1.jacST a1.jacT w1.jacWV a1.jacO3];
dJ = (dy - jac0 * dx)/(dx'*dx) * dx';
figure(7); plot(w0.f,sum(jac1'-jac0'),'r',w0.f,sum(dJ'),'k'); title('ALL')

bt0 = rad2bt(w0.f,w0.rad);  bt1 = rad2bt(w0.f,w1.rad);  dy = bt1-bt0; 
X0 = [pbefore.stemp(iiBin); flipud(pbefore.ptemp(1:97,iiBin)); log(flipud(pbefore.gas_1(1:97,iiBin))); log(flipud(pbefore.gas_2(90,iiBin))); log(flipud(pbefore.gas_3(1:97,iiBin))) ];
X1 = [pafter.stemp(iiBin);  flipud(pafter.ptemp(1:97,iiBin));  log(flipud(pafter.gas_1(1:97,iiBin)));  log(flipud(pafter.gas_2(90,iiBin))); log(flipud(pafter.gas_3(1:97,iiBin))) ];
dx = X1-X0; 
jac0 = [a0.jacST a0.jacT w0.jacWV a0.jacCO2 a0.jacO3];
jac1 = [a1.jacST a1.jacT w1.jacWV a1.jacCO2 a1.jacO3];
dJ = (dy - jac0 * dx)/(dx'*dx) * dx';
figure(7); plot(w0.f,sum(jac1'-jac0'),'r',w0.f,sum(dJ'),'k'); title('ALL')
