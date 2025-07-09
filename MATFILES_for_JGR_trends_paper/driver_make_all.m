% restoredefaultpath       %% testing so there are no hidden paths of mine!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Driver file to loop through figures in JGR paper')
disp('------------------------------------------------');
disp(' ')

for jj = 1 : 4; 
  figure(jj); clf; 
end

disp('Suggest arranging the 4 figures as a 2x2 tile :    Fig 1    |    Fig 2')
disp('                                                   ---------|----------')
disp('                                                   Fig 3    |    Fig 4');
disp('ret to continue'); pause;

for ii = 1 : 14
  clc
  for jj = 1 : 4; 
    figure(jj); clf; 
  end
  str = ['make_jgr_fig' num2str(ii)];
  fprintf(1,'%s \n',str);
  eval(str);
  disp('ret to continue'); pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now do the appendices

clc; 
  for jj = 1 : 4; 
    figure(jj); clf; 
  end
  str = ['figA1_plot'];
  fprintf(1,'%s \n',str);
  eval(str);
  disp('ret to continue'); pause;

clc; 
  for jj = 1 : 4; 
    figure(jj); clf; 
  end
  str = ['figA2_plot'];
  fprintf(1,'%s \n',str);
  eval(str);
  disp('ret to continue'); pause;

clc; 
  for jj = 1 : 4; 
    figure(jj); clf; 
  end
  str = ['make_jgr_figB1'];
  fprintf(1,'%s \n',str);
  eval(str);
  disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
