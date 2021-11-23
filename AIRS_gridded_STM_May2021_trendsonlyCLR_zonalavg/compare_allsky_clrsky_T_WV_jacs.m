load sarta_chans_for_l1c.mat

clrjacKc = zeros(64,2834,300);
cldjacKc = zeros(64,2834,300);

for ii = 1 : 64
  %%%%%%%%%%%%%%%%%%%%%%%%%
  kcrad = ['KCARTAJACS_CLR_T_WV/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = ['load ' kcrad];
  eval(loader)
  rclr(ii,:) = rKc;

  kcjac = ['KCARTAJACS_CLR_T_WV/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat'];
  loader = ['load ' kcjac];
  eval(loader)
  [~,sizeclrjac(ii)] = size(rKc);
    junk = (sizeclrjac(ii)-4)/4;
    indWV = (1:junk)+0*junk; indWVx = (1:junk)+0*100; 
    indO3 = (1:junk)+1*junk; indO3x = (1:junk)+1*100; 
    indTz = (1:junk)+2*junk; indTzx = (1:junk)+2*100; 
  clrjacKc(ii,:,indWVx) = (rKc(:,indWV));
  clrjacKc(ii,:,indO3x) = (rKc(:,indO3));
  clrjacKc(ii,:,indTzx) = (rKc(:,indTz));

  %%%%%%%%%%%%%%%%%%%%%%%%%
  kcrad = ['KCARTAJACS_CLD_T_WV/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = ['load ' kcrad];
  eval(loader)
  rcld(ii,:) = rKc;

  kcjac = ['KCARTAJACS_CLD_T_WV/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat'];
  loader = ['load ' kcjac];
  eval(loader)
  [~,sizecldjac(ii)] = size(rKc);
    junk = (sizecldjac(ii)-4)/6;
    indWV = (1:junk)+0*junk; indWVx = (1:junk)+0*100; 
    indO3 = (1:junk)+1*junk; indO3x = (1:junk)+1*100; 
    indTz = (1:junk)+5*junk; indTzx = (1:junk)+2*100; 
  cldjacKc(ii,:,indWVx) = (rKc(:,indWV));
  cldjacKc(ii,:,indO3x) = (rKc(:,indO3));
  cldjacKc(ii,:,indTzx) = (rKc(:,indTz));
end

subplot(311); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,001:100),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,001:100),1)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('WVjac')
subplot(312); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,101:200),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,101:200),1)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('O3jac')
subplot(313); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,201:300),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,201:300),1)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('Tjac')

subplot(211); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,001:100),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,001:100),1)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('WVjac')
  title('blue : clr    red : allsky');
subplot(212); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,201:300),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,201:300),1)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('Tjac')

subplot(211); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,001:100),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,001:100),1)),2),'r'); plotaxis2; axis([1240 1640 -10 0]); ylabel('WVjac')
  title('blue : clr    red : allsky');
subplot(212); plot(fKc(ichan),sum(squeeze(mean(clrjacKc(:,ichan,201:300),1)),2),'b',fKc(ichan),sum(squeeze(mean(cldjacKc(:,ichan,201:300),1)),2),'r'); plotaxis2; axis([1240 1640 0.75 1]); ylabel('Tjac')

subplot(211); plot(fKc(ichan),sum(squeeze(clrjacKc(32,ichan,001:100)),2),'b',fKc(ichan),sum(squeeze(cldjacKc(32,ichan,001:100)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('WVjac')
  title('blue : clr    red : allsky');
subplot(212); plot(fKc(ichan),sum(squeeze(clrjacKc(32,ichan,201:300)),2),'b',fKc(ichan),sum(squeeze(cldjacKc(32,ichan,201:300)),2),'r'); plotaxis2; xlim([640 1640]); ylabel('Tjac')

subplot(211); plot(fKc(ichan),sum(squeeze(clrjacKc(32,ichan,001:100)),2),'b',fKc(ichan),sum(squeeze(cldjacKc(32,ichan,001:100)),2),'r'); plotaxis2; axis([1240 1640 -10 0]); ylabel('WVjac')
  title('blue : clr    red : allsky');
subplot(212); plot(fKc(ichan),sum(squeeze(clrjacKc(32,ichan,201:300)),2),'b',fKc(ichan),sum(squeeze(cldjacKc(32,ichan,201:300)),2),'r'); plotaxis2; axis([1240 1640 0.75 1]); ylabel('Tjac')

error('sgjsg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
