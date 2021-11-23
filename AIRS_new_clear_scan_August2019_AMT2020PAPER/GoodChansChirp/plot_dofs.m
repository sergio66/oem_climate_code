a_all = load('dofs_cris_airs_chirp.mat');
a_nucaps = load('dofs_cris_airs_chirp_NUCAPS.mat');

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:end-1)-plevs(2:end);
playsD = log(plevs(1:end-1)./plevs(2:end));
plays = playsN./playsD;
plays = plays(4:100);

figure(1);
semilogx(plays,a_all.dofs_T_airsAVG.cdofs,'b',plays,a_all.ndofs_T_airsAVG.cdofs,'bx-',...
         plays,a_nucaps.dofs_T_airsAVG.cdofs,'bd-',plays,a_nucaps.ndofs_T_airsAVG.cdofs,'bs-','linewidth',2)
hold on
semilogx(plays,a_all.dofs_WV_airsAVG.cdofs,'r',plays,a_all.ndofs_WV_airsAVG.cdofs,'rx-',...
         plays,a_nucaps.dofs_WV_airsAVG.cdofs,'rd-',plays,a_nucaps.ndofs_WV_airsAVG.cdofs,'rs-','linewidth',2)
semilogx(plays,a_all.dofs_O3_airsAVG.cdofs,'g',plays,a_all.ndofs_O3_airsAVG.cdofs,'gx-',...
         plays,a_nucaps.dofs_O3_airsAVG.cdofs,'gd-',plays,a_nucaps.ndofs_O3_airsAVG.cdofs,'gs-','linewidth',2)
hold off
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('INSTR allchan','INSTR allchan+noise','INSTR NUCAPS','INSTR NUCAPS+noise','location','best','fontsize',10);
xlabel('P(mb)'); ylabel('DOFS')
title('AIRS')

figure(2);
semilogx(plays,a_all.dofs_T_crisFSR_AVG.cdofs,'b',plays,a_all.ndofs_T_crisFSR_AVG.cdofs,'bx-',...
         plays,a_nucaps.dofs_T_crisFSR_AVG.cdofs,'bd-',plays,a_nucaps.ndofs_T_crisFSR_AVG.cdofs,'bs-','linewidth',2)
hold on
semilogx(plays,a_all.dofs_WV_crisFSR_AVG.cdofs,'r',plays,a_all.ndofs_WV_crisFSR_AVG.cdofs,'rx-',...
         plays,a_nucaps.dofs_WV_crisFSR_AVG.cdofs,'rd-',plays,a_nucaps.ndofs_WV_crisFSR_AVG.cdofs,'rs-','linewidth',2)
semilogx(plays,a_all.dofs_O3_crisFSR_AVG.cdofs,'g',plays,a_all.ndofs_O3_crisFSR_AVG.cdofs,'gx-',...
         plays,a_nucaps.dofs_O3_crisFSR_AVG.cdofs,'gd-',plays,a_nucaps.ndofs_O3_crisFSR_AVG.cdofs,'gs-','linewidth',2)
hold off
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('INSTR allchan','INSTR allchan+noise','INSTR NUCAPS','INSTR NUCAPS+noise','location','best','fontsize',10);
xlabel('P(mb)'); ylabel('DOFS')
title('CRIS FSR')

figure(3);
semilogx(plays,a_all.dofs_T_chirpA_AVG.cdofs,'b',plays,a_all.ndofs_T_chirpA_AVG.cdofs,'bx-',...
         plays,a_nucaps.dofs_T_chirpA_AVG.cdofs,'bd-',plays,a_nucaps.ndofs_T_chirpA_AVG.cdofs,'bs-','linewidth',2)
hold on
semilogx(plays,a_all.dofs_WV_chirpA_AVG.cdofs,'r',plays,a_all.ndofs_WV_chirpA_AVG.cdofs,'rx-',...
         plays,a_nucaps.dofs_WV_chirpA_AVG.cdofs,'rd-',plays,a_nucaps.ndofs_WV_chirpA_AVG.cdofs,'rs-','linewidth',2)
semilogx(plays,a_all.dofs_O3_chirpA_AVG.cdofs,'g',plays,a_all.ndofs_O3_chirpA_AVG.cdofs,'gx-',...
         plays,a_nucaps.dofs_O3_chirpA_AVG.cdofs,'gd-',plays,a_nucaps.ndofs_O3_chirpA_AVG.cdofs,'gs-','linewidth',2)
hold off
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('INSTR allchan','INSTR allchan+noise','INSTR NUCAPS','INSTR NUCAPS+noise','location','best','fontsize',10);
xlabel('P(mb)'); ylabel('DOFS')
title('CHIRP A')

figure(4);
semilogx(plays,a_all.dofs_T_chirpC_AVG.cdofs,'b',plays,a_all.ndofs_T_chirpC_AVG.cdofs,'bx-',...
         plays,a_nucaps.dofs_T_chirpC_AVG.cdofs,'bd-',plays,a_nucaps.ndofs_T_chirpC_AVG.cdofs,'bs-','linewidth',2)
hold on
semilogx(plays,a_all.dofs_WV_chirpC_AVG.cdofs,'r',plays,a_all.ndofs_WV_chirpC_AVG.cdofs,'rx-',...
         plays,a_nucaps.dofs_WV_chirpC_AVG.cdofs,'rd-',plays,a_nucaps.ndofs_WV_chirpC_AVG.cdofs,'rs-','linewidth',2)
semilogx(plays,a_all.dofs_O3_chirpC_AVG.cdofs,'g',plays,a_all.ndofs_O3_chirpC_AVG.cdofs,'gx-',...
         plays,a_nucaps.dofs_O3_chirpC_AVG.cdofs,'gd-',plays,a_nucaps.ndofs_O3_chirpC_AVG.cdofs,'gs-','linewidth',2)
hold off
set(gca,'xdir','reverse'); grid
xlim([1 1000])
hl = legend('INSTR allchan','INSTR allchan+noise','INSTR NUCAPS','INSTR NUCAPS+noise','location','best','fontsize',10);
xlabel('P(mb)'); ylabel('DOFS')
title('CHIRP C')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
