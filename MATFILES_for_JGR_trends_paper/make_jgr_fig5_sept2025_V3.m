figfileA = '/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5a_L1Cspectraltrends_bigfont.fig';
figfileB = '/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5b_L1Cspectraltrends_unc_bigfont.fig';
figfileC = '/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5c_ERA5spectraltrends_bigfont.fig';

figure(1); clf; % close
  loader = ['hgload ' figfileA ';']; eval(loader) 
figure(2); clf; % close
  loader = ['hgload ' figfileB ';']; eval(loader) 
figure(3); clf; % close
  loader = ['hgload ' figfileC ';']; eval(loader) 


