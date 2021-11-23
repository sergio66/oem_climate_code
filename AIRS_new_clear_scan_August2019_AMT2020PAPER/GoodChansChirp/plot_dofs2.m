a_all = load('dofs_cris_airs_chirp.mat');
a_nucaps = load('dofs_cris_airs_chirp_NUCAPS.mat');

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:end-1)-plevs(2:end);
playsD = log(plevs(1:end-1)./plevs(2:end));
plays = playsN./playsD;
plays = plays(4:100);

figure(1);
semilogx(plays,a_all.ndofs_T_airsAVG.cdofs,'b',plays,a_all.ndofs_T_crisFSR_AVG.cdofs,'r',...
         plays,a_all.ndofs_T_chirpA_AVG.cdofs,'k',plays,a_all.ndofs_T_chirpC_AVG.cdofs,'g',...
         plays,a_nucaps.ndofs_T_airsAVG.cdofs,'bo-',plays,a_nucaps.ndofs_T_crisFSR_AVG.cdofs,'ro-',...
         plays,a_nucaps.ndofs_T_chirpA_AVG.cdofs,'ko-',plays,a_nucaps.ndofs_T_chirpC_AVG.cdofs,'go-',...
         'linewidth',2)
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('AIRS allchan','CRIS FSR allchan','CHIRP A allchan','CHIRP C allchan',...
            'AIRS NUCAPS','CRIS FSR NUCAPS','CHIRP A NUCAPS','CHIRP C NUCAPS',...
            'location','best','fontsize',10);
xlabel('P(mb)'); ylabel('T DOFS (0.2K noise)')
title('Tz')

figure(2);
semilogx(plays,a_all.ndofs_WV_airsAVG.cdofs,'b',plays,a_all.ndofs_WV_crisFSR_AVG.cdofs,'r',...
         plays,a_all.ndofs_WV_chirpA_AVG.cdofs,'k',plays,a_all.ndofs_WV_chirpC_AVG.cdofs,'g',...
         plays,a_nucaps.ndofs_WV_airsAVG.cdofs,'bo-',plays,a_nucaps.ndofs_WV_crisFSR_AVG.cdofs,'ro-',...
         plays,a_nucaps.ndofs_WV_chirpA_AVG.cdofs,'ko-',plays,a_nucaps.ndofs_WV_chirpC_AVG.cdofs,'go-',...
         'linewidth',2)
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('AIRS allchan','CRIS FSR allchan','CHIRP A allchan','CHIRP C allchan',...
            'AIRS NUCAPS','CRIS FSR NUCAPS','CHIRP A NUCAPS','CHIRP C NUCAPS',...
            'location','best','fontsize',10);
xlabel('P(mb)'); ylabel('WV DOFS (0.2K noise)')
title('WV')

figure(3);
semilogx(plays,a_all.ndofs_O3_airsAVG.cdofs,'b',plays,a_all.ndofs_O3_crisFSR_AVG.cdofs,'r',...
         plays,a_all.ndofs_O3_chirpA_AVG.cdofs,'k',plays,a_all.ndofs_O3_chirpC_AVG.cdofs,'g',...
         plays,a_nucaps.ndofs_O3_airsAVG.cdofs,'bo-',plays,a_nucaps.ndofs_O3_crisFSR_AVG.cdofs,'ro-',...
         plays,a_nucaps.ndofs_O3_chirpA_AVG.cdofs,'ko-',plays,a_nucaps.ndofs_O3_chirpC_AVG.cdofs,'go-',...
         'linewidth',2)
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('AIRS allchan','CRIS FSR allchan','CHIRP A allchan','CHIRP C allchan',...
            'AIRS NUCAPS','CRIS FSR NUCAPS','CHIRP A NUCAPS','CHIRP C NUCAPS',...
            'location','best','fontsize',10);
xlabel('P(mb)'); ylabel('O3 DOFS (0.2K noise)')
title('O3')

for ii = 1 : 3
  figure(ii); grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPrint = -1;
iPrint = +1;
if iPrint > 0
  addpath /asl/matlib/plotutils
  disp('saving fig files')
  figure(1); aslprint('dofsTz.pdf');
  figure(2); aslprint('dofsWV.pdf');
  figure(3); aslprint('dofsO3.pdf');
end
