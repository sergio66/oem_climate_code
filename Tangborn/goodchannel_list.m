goody = good;
if (good == -1) 
  g = 1 : 2378; 
elseif (good == 1)          %%june 14 
  good = '/home/strow/Matlab/Airs/igood_jun14.mat';  
  loader = ['load  ' good]; 
  eval(loader); 
  g = igood; 
elseif (good == 2)          %%july 20 
  good = '/home/strow/Matlab/Airs/clist_v6.0.4.mat'; 
  loader = ['load  ' good]; 
  eval(loader); 
  g = find(bad == 0); 
elseif (good == 3)          %%aug 31 
  good = '/home/strow/Matlab/Airs/clist_v6.0.4.mat'; 
  good = '/asl/data/airs/srf/clist.mat'; 
  loader = ['load  ' good]; 
  eval(loader); 
  bad(1791) = 1; 
  g = find(bad == 0); 
  
  %%use ibad.mat upto Jan 15, 2003
  if (year==2002)
    clear g
%    good = '/yam/s1/asl/ops/airs/tlscf/2002/10/ibad.mat';
    good = 'ibad.mat';
    loader = ['load  ' good]; 
    eval(loader); 
    g = 1:2378;
    g = setdiff(g,ibad);

%  good = '/home/strow/Matlab/Airs/clist_v6.0.4.mat'; 
  good = '/asl/data/airs/srf/clist.mat'; 
  loader = ['load  ' good]; 
  eval(loader); 
  bad(1791) = 1; 
  g = find(bad == 0); 

  elseif ((year==2003) & (month == 1) & (date <= 15))
    clear g
%    good = '/yam/s1/asl/ops/airs/tlscf/2002/10/ibad.mat';
    good = 'ibad.mat';
    loader = ['load  ' good]; 
    eval(loader); 
    g = 1:2378;
    g = setdiff(g,ibad);

  %%use ibad2.mat from  Jan 15, 2003 to Oct 31, 2003
  elseif ((year==2003) & (month == 1) & (date >= 15))
    clear g
%    good = '/yam/s1/asl/ops/airs/tlscf/2003/03/ibad2.mat';
    good = 'ibad2.mat';
    loader = ['load  ' good]; 
    eval(loader); 
    g = 1:2378;
    g = setdiff(g,ibad2); 
  elseif ((year==2003) & (month <= 10))
    clear g
%    good = '/yam/s1/asl/ops/airs/tlscf/2003/03/ibad2.mat';
    good = 'ibad2.mat';
    loader = ['load  ' good]; 
    eval(loader); 
    g = 1:2378;
    g = setdiff(g,ibad2);
    end
  g=g';
else 
  error('invalid "good" parameter') 
  end 

%good = '/asl/data/airs/srf/master_ig.mat';
good = 'master_ig.mat';
loader = ['load  ' good]; 
eval(loader); 
g = ig;
g = g';

goodchannels = g;
good = goody;
