function [] = hpc_usage(Q,ocb_set,iLinearOrLog)

% hpc_usage(Q,ocb_set,iLinearOrLog)
%  defaults Q = 16, ocb_set = 0,iLinearOrLog = 1

save_no_legend_autoupdate

disp('checking hpc nodes health ...')

if nargin == 0
  Q = 16;
  ocb_set = 0;
  iLinearOrLog = 1;
elseif nargin == 1
  ocb_set = 0;
  iLinearOrLog = 1;
elseif nargin == 2
  iLinearOrLog = 1;
end

%% alias sqrun='squeue -u sergio -t R -o "%.19i %.9P %.8j %.8u %.2t %.10M %.6D %R"'
sqrun = 'squeue -u sergio -t R -o "%.19i %.9P %.8j %.8u %.2t %.10M %.6D %R"';

if ocb_set == 0
  thedir = dir(['Output/Quantile' num2str(Q,'%02d') '/test*.mat']);
  fprintf(1,'looking in Output/Quantile %2i /test*.mat \n',Q);
elseif ocb_set == 1
  thedir = dir(['Output_CAL/Quantile' num2str(Q,'%02d') '/test*.mat']);
  fprintf(1,'looking in Output_CAL/Quantile %2i /test*.mat \n',Q);
end

eval(['!' sqrun ' >& /home/sergio/ugh']);

%{
expecting
0        1         2         3         4         5         6         7         8         9        10        11        12       13         14
1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
         11443649_1  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_2  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_3  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_4  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_5  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_6  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_7  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_8  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
         11443649_9  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
        11443649_10  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
        11443649_11  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
        11443649_12  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
        11443649_13  high_mem  ANOMALY   sergio  R       1:41      1 cnode013
%}

iarray = 0;
inode = 0;
jobsPERproessor = 72; 

ii   = 0;
fid = fopen('/home/sergio/ugh','r');
tline = fgetl(fid);
while ischar(tline)
  ii = ii + 1;
  if ii >= 2
    xline{ii-1} = tline;
    
    xpart_name  = tline(22:30);
    xarrayID   = tline(7:19);
      boo = findstr(xarrayID,'_');
      xarrayID = str2num(xarrayID(boo+1:end));
    xstatus     = tline(50:50);
    xnode_name  = tline(70:end);
      xnode_name = str2num(xnode_name(6:end));
    %fprintf(1,'arrayID %4i cnode%4i \n',xarrayID,xnode_name);

    iS = (xarrayID-1)*jobsPERproessor+1;
    iE = (xarrayID-1)*jobsPERproessor+jobsPERproessor;
    
    iFound = 0;
    for kk = iS : iE
      if ocb_set == 0
        fname = ['Output/Quantile' num2str(Q,'%02d') '/test' num2str(kk) '.mat'];
      else
        fname = ['Output_CAL/Quantile' num2str(Q,'%02d') '/test' num2str(kk) '.mat'];
      end
      if exist(fname)
         iFound = iFound + 1;
       end
     end
    fprintf(1,' cnode%3i jobarray %3i %04i-%04i found %2i \n',xnode_name,xarrayID,iS,iE,iFound);

    iarray = xarrayID;
    node(iarray,1) = xnode_name;
    node(iarray,2) = xarrayID;
    node(iarray,3) = iFound;

    if inode == 0
      inode = inode + 1;
      expected(inode,1) = xnode_name;
      expected(inode,2) = iFound;
      expected(inode,3) = jobsPERproessor;
    elseif inode > 0
      [moo,iA,iB] = intersect(expected(:,1),xnode_name);
      if length(moo) == 0
        %% new node
        inode = inode + 1;
        expected(inode,1) = xnode_name;
        expected(inode,2) = iFound;
        expected(inode,3) = jobsPERproessor;
      elseif length(moo) == 1
        %% node exists
        expected(iA,1) = xnode_name;
        expected(iA,2) = iFound + expected(iA,2);
        expected(iA,3) = jobsPERproessor + expected(iA,3);
      end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end  %% if ii > 1

  tline = fgetl(fid);
end

if exist('iFound')
  figure(1); plot(node(:,1),node(:,2),'o'); xlabel('cnode#'); title('jobs done');
  Yi = find(node(:,1) > 0);  %% if nodes eg cnode031 and cnode054 are the only one running, then it fills in stuff with ZEROS
  jettY = unique(node(Yi,1));
  jettY = min(jettY) : max(jettY);
  if length(jettY) == 1
    jettY = jettY-1:jettY+1;
  end
  jettY = colormap(jet(length(jettY)));
  figure(1); scatter(node(Yi,2),node(Yi,3),40,node(Yi,1),'filled'); colorbar; colormap(jettY);
    xlabel('array ID (1-64)'); ylabel(['jobs done out of ' num2str(jobsPERproessor)]); colorbar; title('colorbar = cnode')
  
  [Y,I] = sort(expected(:,1));
  figure(2); plot(expected(I,1),expected(I,2) ./ expected(I,3),'-o');
    xlabel('node number'); ylabel(['fraction jobs done (out of ' num2str(jobsPERproessor) ')']);
    grid
  wah = [expected(I,1)    expected(I,2) ./ expected(I,3)];
  fprintf(1,'cnode %3i   fraction done %8.4f \n',wah')
else
  disp(' oooer nothing running')
end

if ocb_set == 0
  fprintf(1,'found  %4i files in Output/Quantile%2i \n',length(thedir),Q);
elseif ocb_set == 1
  fprintf(1,'found  %4i files in Output_CAL/Quantile%2i \n',length(thedir),Q);
end
