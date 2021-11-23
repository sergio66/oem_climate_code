%% this is to model radiances vs jacs
%% works quite well for low OD, pretty crappy for large OD

mu = cos(p.scanang*pi/180);

%% quasi window channel eg 780 cm-1
f0 = 780;
T0 = 300; T1 = 260;           %% surf temp, and temp of emitting layer
tau_all=1; tau_emit = 0.25;   %% tau_all = ground to space OD, tau_emit = peak of wgt_fcn to space OD

%% window channel eg 900 cm-1
f0 = 900;
T0 = 300; T1 = 290;           %% surf temp, and temp of emitting layer
tau_all=1; tau_emit = 0.1;    %% tau_all = ground to space OD, tau_emit = peak of wgt_fcn to space OD

%% high alt sounding channel eg 710 cm-1
f0 = 710;
T0 = 300; T1 = 230;           %% surf temp, and temp of emitting layer
tau_all=10; tau_emit = 1;     %% tau_all = ground to space OD, tau_emit = peak of wgt_fcn to space OD

[T0 T1]

figure(1); plot(p.scanang,1./mu,'b.-',p.scanang,exp(-tau_all./mu),'r.-'); 
hl = legend('jac(\theta)','rad(\theta)'); xlabel('scanang'); title('Very Very Simple Model')

figure(2); plot(p.scanang,ttorad(f0,T0)*exp(-tau_all./mu)*1./mu,'b.-',p.scanang,ttorad(f0,T0)*exp(-tau_all./mu),'r.-'); 
hl = legend('jac(\theta)','rad(\theta)'); xlabel('scanang'); title('Very Simple Model')

radx = ttorad(f0,T0)*exp(-tau_all./mu) + (1-exp(-tau_emit./mu)).*exp(-tau_emit./mu)*ttorad(f0,T1);
jacx = -ttorad(f0,T0)*1./mu .*exp(-tau_all./mu) + -1./mu.*(1-exp(-tau_emit./mu)).*1./mu.*exp(-tau_emit./mu)*ttorad(f0,T1);
figure(3); plot(p.scanang,-jacx,'b.-',p.scanang,radx,'r.-'); hl = legend('jac(\theta)','rad(\theta)'); xlabel('scanang'); title('Better Model')
