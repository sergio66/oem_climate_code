load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q01.mat
moo01 = b_desc;
load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q02.mat
moo02 = b_desc;
load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat
moo03 = b_desc;
load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q04.mat
moo04 = b_desc;
load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat
moo05 = b_desc;
load h2645structure.mat

zoo01 = permute(moo01,[3 1 2]); zoo01 = reshape(zoo01,2645,72*64);
zoo02 = permute(moo02,[3 1 2]); zoo02 = reshape(zoo02,2645,72*64);
zoo03 = permute(moo03,[3 1 2]); zoo03 = reshape(zoo03,2645,72*64);
zoo04 = permute(moo04,[3 1 2]); zoo04 = reshape(zoo04,2645,72*64);
zoo05 = permute(moo05,[3 1 2]); zoo05 = reshape(zoo05,2645,72*64);

figure(1)
plot(1:2645,nanmean(zoo01,2))
plot(h.vchan,nanmean(zoo01,2),h.vchan,nanmean(zoo02,2))
plot(h.vchan,nanmean(zoo01,2),h.vchan,nanmean(zoo02,2),h.vchan,nanmean(zoo03,2),h.vchan,nanmean(zoo04,2),h.vchan,nanmean(zoo05,2))
grid;
xlim([640 1640])
hl = legend('Q1','Q2','Q3','Q4','Q5','location','best','fontsize',8);

woof = nanmean(zoo01,2);
subplot(211); plot(h.vchan,woof); ylabel('Q01'); grid; xlim([640 1640])
subplot(212); plot(h.vchan,woof-nanmean(zoo01,2),h.vchan,woof-nanmean(zoo02,2),h.vchan,woof-nanmean(zoo03,2),h.vchan,woof-nanmean(zoo04,2),h.vchan,woof-nanmean(zoo05,2)); ylabel('Q01-Qx')
grid;
xlim([640 1640]); ax = axis; ylim([-1e-3 ax(4)])
hl = legend('Q1','Q2','Q3','Q4','Q5','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

