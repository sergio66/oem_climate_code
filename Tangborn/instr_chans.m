function y = instr_chans(instr,iFreqOrNoise)

% iFreqOrNoise = 1 for freq, 2 for NEdT

if nargin == 0
  instr = 'airs';
  iFreqOrNoise = 1;  %% gets freqs
end

if nargin == 1
  iFreqOrNoise = 1;  %% gets freqs
end

if iFreqOrNoise == 1
  if findstr(instr,'iasi') | findstr(instr,'IASI')
    y = 1:8461;
    y = 645 + (y-1)*0.25;
  elseif findstr(instr,'cris') | findstr(instr,'CRIS')
    fo = [(650-2*0.625):0.625:(1095+2*0.625),  ...
          (1210-2*1.25):1.25:(1750+2*1.25), ...
          (2155-2*2.5):2.5:(2550+2*2.5)]';
    load cris_chans_jan2012.mat
    y = hcris.vchan;
  elseif findstr(instr,'airs') | findstr(instr,'AIRS')
    fx = 'airs_chanlist_calib_tomp';

    dd = load(fx);

    %% everything is screwup here ie need to re-order things 
    %% from this Lockhed Martin file
    dd2 = dd(:,2); [Y,I] = sort(dd2);
    dd = dd(I,:);
    y = dd(:,3); 
  end
end

if iFreqOrNoise == 2
  if findstr(instr,'iasi') | findstr(instr,'IASI')
    y = 1:8461;
    y = 645 + (y-1)*0.25;
    error('iasi noise???')
  elseif findstr(instr,'cris') | findstr(instr,'CRIS')
    fo = [(650-2*0.625):0.625:(1095+2*0.625),  ...
          (1210-2*1.25):1.25:(1750+2*1.25), ...
          (2155-2*2.5):2.5:(2550+2*2.5)]';
    load cris_chans_jan2012.mat
    y = hcris.vchan;
    error('cris noise???')
  elseif findstr(instr,'airs') | findstr(instr,'AIRS')
    fx = 'AIRS_IASI_AMSU_JACS_gas_T/airs_chanlist_calib_tomp';
    dd = load(fx);
    %% everything is screwup here ie need to re-order things 
    %% from this Lockhed Martin file
    dd2 = dd(:,2); [Y,I] = sort(dd2);
    dd = dd(I,:);
    y = dd(:,6); 
  end
end
