function new_m_ts_jac0 = combine_sarta_anaytic_precomputed_jacs(xm_ts_jac0,xnlays,xqrenorm,m_ts_jac0,nlays,qrenorm);

[xmm,xnn] = size(xm_ts_jac0);
[ mm, nn] = size(m_ts_jac0);

  %% /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/LoopMake_CFC_JacskCARTA/compare_sarta_kcarta_jac.m

%% do the CO2/N2O/CH4 jacs
xm_ts_jac0(:,1) = xm_ts_jac0(:,1)*2.2/370;
xm_ts_jac0(:,2) = xm_ts_jac0(:,2)*1/300;
xm_ts_jac0(:,3) = xm_ts_jac0(:,3)*5/1860;

%% do the CFC jacs
%xm_ts_jac0(:,4) = m_ts_jac0(:,4) * 1/qrenorm(4);
%xm_ts_jac0(:,5) = m_ts_jac0(:,5) * 1/qrenorm(5);
xm_ts_jac0(:,4) = m_ts_jac0(:,4);
xm_ts_jac0(:,5) = m_ts_jac0(:,5);

%% do the ST jac
xm_ts_jac0(:,6) = xm_ts_jac0(:,6) * qrenorm(6);

if xnn ~= 6+xnlays*3
  [xnn xnlays]
  error('combine_sarta_anaytic_precomputed_jacs.m : xnlays ~= xnn')
end
if nn ~= 6+nlays*3
  [nn nlays]
  error('combine_sarta_anaytic_precomputed_jacs.m : nlays ~= nn')
end

if xnlays ~= nlays & xnn == nn
  error('combine_sarta_anaytic_precomputed_jacs.m : xnlays ~= nlays & xnn == nn')
elseif xnlays == nlays & xnn ~= nn
  error('combine_sarta_anaytic_precomputed_jacs.m : xnlays == nlays & xnn ~= nn')

elseif xnlays == nlays & xnn == nn
  numlays = (nn-6)/3;
  ind = (1:numlays) + 6;  xm_ts_jac0(:,ind) = xm_ts_jac0(:,ind) * mean(qrenorm(ind));  %% WV
  ind = ind + numlays;    xm_ts_jac0(:,ind) = xm_ts_jac0(:,ind) * mean(qrenorm(ind));  %% T
  ind = ind + numlays;    xm_ts_jac0(:,ind) = xm_ts_jac0(:,ind) * mean(qrenorm(ind));  %% O3

elseif xnlays > nlays

  numlays  = (nn-6)/3;
  xnumlays = (xnn-6)/3;

  junk = [];
  junk(:,1:6) = xm_ts_jac0(:,1:6);
  ind = (1:numlays) + 6;  xindx = (1:xnumlays) + 6;  xind = xindx(1:numlays); junk(:,ind) = xm_ts_jac0(:,xind) * mean(qrenorm(ind));  %% WV
  ind = ind + numlays;    xindx = xindx + xnumlays;  xind = xindx(1:numlays); junk(:,ind) = xm_ts_jac0(:,xind) * mean(qrenorm(ind));  %% T
  ind = ind + numlays;    xindx = xindx + xnumlays;  xind = xindx(1:numlays); junk(:,ind) = xm_ts_jac0(:,xind) * mean(qrenorm(ind));  %% O3

  fprintf(1,'size(kcarta jac we want) = %5i x %5i \n',mm,nn);
  fprintf(1,'size(sarta analytic jac) = %5i x %5i \n',xmm,xnn);
  [xmmx,xnnx] = size(junk);
  fprintf(1,'size(new   analytic jac) = %5i x %5i \n',xmmx,xnnx);
  xm_ts_jac0 = junk;

elseif xnlays < nlays

  numlays  = (nn-6)/3;
  xnumlays = (xnn-6)/3;

  dn = abs(numlays - xnumlays);

  junk = [];
  junk(:,1:6) = xm_ts_jac0(:,1:6);
  ind = (1:numlays) + 6;  xindx = (1:xnumlays) + 6;  boo = ind(1:xnumlays); junk(:,boo) = xm_ts_jac0(:,xindx); junk(:,boo(end)+1:ind(end)) = junk(:,boo(end)) * ones(1,dn); junk(:,ind) = junk(:,ind) * mean(qrenorm(ind));  %% WV
  ind = ind + numlays;    xindx = xindx + xnumlays;  boo = ind(1:xnumlays); junk(:,boo) = xm_ts_jac0(:,xindx); junk(:,boo(end)+1:ind(end)) = junk(:,boo(end)) * ones(1,dn); junk(:,ind) = junk(:,ind) * mean(qrenorm(ind));  %% T
  ind = ind + numlays;    xindx = xindx + xnumlays;  boo = ind(1:xnumlays); junk(:,boo) = xm_ts_jac0(:,xindx); junk(:,boo(end)+1:ind(end)) = junk(:,boo(end)) * ones(1,dn); junk(:,ind) = junk(:,ind) * mean(qrenorm(ind));  %% O3

  fprintf(1,'size(kcarta jac we want) = %5i x %5i \n',mm,nn);
  fprintf(1,'size(sarta analytic jac) = %5i x %5i \n',xmm,xnn);
  [xmmx,xnnx] = size(junk);
  fprintf(1,'size(new   analytic jac) = %5i x %5i \n',xmmx,xnnx);
  xm_ts_jac0 = junk;
end

bad = find(isnan(xm_ts_jac0) | isinf(xm_ts_jac0));
xm_ts_jac0(bad) = 0;
new_m_ts_jac0 = xm_ts_jac0;

iDebugPlot = +1;
iDebugPlot = -1;
if iDebugPlot > 0
  figure(4); clf; plot(1:2645,new_m_ts_jac0(:,6),1:2645,m_ts_jac0(:,6)); legend('new','old','location','best','fontsize',10); title('Stemp')
  figure(5); clf; plot(1:2645,new_m_ts_jac0(:,1),1:2645,m_ts_jac0(:,1)); legend('new','old','location','best','fontsize',10); title('CO2')
  figure(6); clf; plot(1:2645,new_m_ts_jac0(:,2),1:2645,m_ts_jac0(:,2)); legend('new','old','location','best','fontsize',10); title('N2O')
  figure(7); clf; plot(1:2645,new_m_ts_jac0(:,3),1:2645,m_ts_jac0(:,3)); legend('new','old','location','best','fontsize',10); title('CH4')
  %plot(1:2645,new_m_ts_jac0(:,4),1:2645,m_ts_jac0(:,4))
  %plot(1:2645,new_m_ts_jac0(:,4),'.',1:2645,m_ts_jac0(:,4))
  %plot(1:2645,new_m_ts_jac0(:,5),'.',1:2645,m_ts_jac0(:,5))
  %plot(1:2645,new_m_ts_jac0(:,6),1:2645,m_ts_jac0(:,6))

  figure(1); clf; ind = (1:97)+6 + 0*97; plot(1:2645,sum(new_m_ts_jac0(:,ind),2),1:2645,sum(m_ts_jac0(:,ind),2)); legend('new','old','location','best','fontsize',10); title('col WV')
  figure(2); clf; ind = (1:97)+6 + 1*97; plot(1:2645,sum(new_m_ts_jac0(:,ind),2),1:2645,sum(m_ts_jac0(:,ind),2)); legend('new','old','location','best','fontsize',10); title('col TZ')
  figure(3); clf; ind = (1:97)+6 + 2*97; plot(1:2645,sum(new_m_ts_jac0(:,ind),2),1:2645,sum(m_ts_jac0(:,ind),2)); legend('new','old','location','best','fontsize',10); title('col O3')

end

