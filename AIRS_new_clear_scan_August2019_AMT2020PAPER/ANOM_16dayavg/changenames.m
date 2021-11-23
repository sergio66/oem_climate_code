for ii = 1 : 40
  mver = ['!mv latbin' num2str(ii) '.mat latbin_180dayavg_' num2str(ii) '.mat'];
  eval(mver)

  mver = ['!mv latbin' num2str(ii) '_cal.mat latbin_180dayavg_' num2str(ii) '_cal.mat'];
  eval(mver)
end
