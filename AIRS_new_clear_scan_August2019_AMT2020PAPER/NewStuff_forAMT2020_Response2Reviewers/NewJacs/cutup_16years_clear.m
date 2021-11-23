addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools/
addpath /asl/matlib/rtptools/
addpath /asl/matlib/gribtools/

iTotal = 0;
for yy = 2002 : 2019
  dirx = ['/asl/rtp/rtp_airicrad_v6/clear/' num2str(yy,'%4d')];
  daysINmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
  if mod(yy,4) == 0
    daysINmonth(2) = 29;
  end
  mmStart = 1; mmEnd = 12;
  if yy == 2002
    mmStart = 9;
  elseif yy == 2019
    mmStart = 8;
  end

  for mm = mmStart : mmEnd
    iFound = -1;
    iCnt = 0;
    while iFound < 0 & iCnt < 30
      randday = floor(30*rand(1))+1;
      iCnt = iCnt + 1;
      if mm > 1
        doy = sum(daysINmonth(1:mm-1)) + randday;
      else
        doy = randday;
      end
      fname = [dirx '/era_airicrad_day' num2str(doy,'%03d') '_clear.rtp'];
      if exist(fname)
        iFound = +1;
      end
    end    

    iTotal = iTotal + 1;
    fprintf(1,'%3i : %4i %2i %s \n',iTotal,yy,mm,fname);
    [h,ha,p,pa] = rtpread(fname);
    if ~isfield(p,'rcalc') & isfield(p,'rclr')
      p.rcalc = p.rclr;
      p = rmfield(p,'rclr');
    end

    fprintf(1,'  have %5i fovs in here \n',length(p.stemp))

    boo = find(p.solzen > 90 & p.landfrac == 0);
    [h,p] = subset_rtp(h,p,[],[],boo);
    fprintf(1,'  of which %5i fovs are night/ocean \n',length(p.stemp))

    %moo = unique(floor(750*rand(1,min(length(p.stemp),750)))+1);
    moo = rand(1,750)*(length(p.stemp)-1);
    moo = floor(moo) + 1;
    moo = unique(moo);

    [h,p] = subset_rtp(h,p,[],[],moo);
    fprintf(1,'  of which we select %5i fovs \n',length(p.stemp))

    fout = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2002/08/31/cutup_16years_clear_' num2str(yy,'%04d') '_' num2str(mm,'%02d') '_' num2str(iTotal,'%03d') '.rtp'];
    fout = ['/asl/s1/sergio/rtp/rtp_airicrad_v6/2002/08/31/cutup_16years_clear_' num2str(iTotal,'%03d') '.rtp'];
    rtpwrite(fout,h,ha,p,pa);
  end
end

