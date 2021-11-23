addpath /home/sergio/KCARTA/MATLAB

iaFound = zeros(1,40);

iWhich = -1; %% before
iWhich = +1; %% after

for ii = 1 : 40
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  radIN = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  jacIN = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_coljac.mat'];
  if exist(jacIN) & exist(radIN)
    iaFound(ii) = 1;

    loader = ['a = load(''' radIN ''');'];
    eval(loader);
    rad = a.rKc;

    loader = ['a = load(''' jacIN ''');'];
    eval(loader);
    jac = a.rKc;

    ind = 1; jacCO2   = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 2; jacN2O   = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 3; jacCO    = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 4; jacCH4   = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 5; jacCFC11 = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 6; jacCFC12 = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 7; jacT     = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    ind = 8; jacST    = rad2bt(a.fKc,jac(:,ind))-rad2bt(a.fKc,rad); 
    clear jac
    f = a.fKc;
    plot(f,jacCO2,'b',f,jacN2O,'g',f,jacCH4,'r',f,jacST,'k'); pause(0.1);

    if iWhich == +1
      fout = ['PAFTER/coljacresults' num2str(ii,'%02d') '.mat'];
    elseif iWhich == -1
      fout = ['PBEFORE/coljacresults' num2str(ii,'%02d') '.mat'];
    end

    if exist(fout)
      fprintf(1,'%s already exists ... not saving \n',fout);
    else
      saver = ['save ' fout ' f rad jac*'];
      eval(saver)
    end
  
    pause(0.1)
  end
end

fprintf(1,'\n');
 
if length(find(iaFound == 0)) > 0
  disp('oops not done all!')
  find(iaFound == 0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
a0 = load('PBEFORE/jacresults01_G_2356.mat');
c0 = load('PBEFORE/coljacresults01.mat');
a1 = load('PAFTER/jacresults01_G_2356.mat');
c1 = load('PAFTER/coljacresults01.mat');
plot(a0.f,a0.jacCO2/100,'b.',c0.f,c0.jacCO2/0.1,'r')
plot(a0.f,a1.jacCO2/100 - a0.jacCO2/100,'b.',c0.f,c1.jacCO2/0.1 - c0.jacCO2/0.1,'r')
%}
