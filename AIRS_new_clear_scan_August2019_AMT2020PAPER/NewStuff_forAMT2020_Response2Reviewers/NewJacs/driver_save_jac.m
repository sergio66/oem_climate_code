addpath /home/sergio/KCARTA/MATLAB

iaFound = zeros(1,40);

iWhich = +1; %% after
iWhich = -1; %% before

iG = 2356;   %% 2356
iG = 1;      %% WV

for ii = 1 : 40
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  radIN = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  jacIN = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_jac.mat'];
  if exist(jacIN) & exist(radIN)
    iaFound(ii) = 1;

    loader = ['a = load(''' radIN ''');'];
    eval(loader);
    rad = a.rKc;

    loader = ['a = load(''' jacIN ''');'];
    eval(loader);
    jac = a.rKc;

    if iG == 2356
      ind = (1:97) + 0*97; jacCO2 = jac(:,ind); jacCO2 = sum(jacCO2,2);
      ind = (1:97) + 1*97; jacO3 = jac(:,ind);
      ind = (1:97) + 2*97; jacN2O = jac(:,ind); jacN2O = sum(jacN2O,2);
      ind = (1:97) + 3*97; jacCH4 = jac(:,ind); jacCH4 = sum(jacCH4,2);
      ind = (1:97) + 4*97; jacT   = jac(:,ind); 
      ind = (1:97) + 5*97; wgt    = jac(:,ind);
      ind = max(ind)+1;    jacST  = jac(:,ind);
      clear jac
      f = a.fKc;
      plot(f,jacCO2,'b',f,jacN2O,'g',f,jacCH4,'r',f,jacST,'k'); pause(0.1);
    elseif iG == 1
      ind = (1:97) + 0*97; jacWV = jac(:,ind); 
      ind = (1:97) + 1*97; jacT   = jac(:,ind); 
      ind = (1:97) + 2*97; wgt    = jac(:,ind);
      ind = max(ind)+1;    jacST  = jac(:,ind);
      clear jac
      f = a.fKc;
    end

    if iWhich == +1
      fout = ['PAFTER/jacresults' num2str(ii,'%02d') '_G_' num2str(iG) '.mat'];
    elseif iWhich == -1
      fout = ['PBEFORE/jacresults' num2str(ii,'%02d') '_G_' num2str(iG) '.mat'];
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
