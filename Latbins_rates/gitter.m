clc
str = ['!git add override_defaults_latbins_*.m']; eval(str)
str = ['!git status'];                            eval(str)
lala = clock; 
  datex = ['Updates : ' date ' ' num2str(lala(:,4),'%02d') ':' num2str(lala(:,5),'%02d')];
  fprintf(1,'%s \n',datex);
  str = ['!git commit -m "' datex '" '];
  eval(str)
str = ['!git push origin']; eval(str)
