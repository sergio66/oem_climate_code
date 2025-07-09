%load figa1_data;   % mei_time mei2 tile_time tile_anomaly
load figA1_data;    % mei_time mei2 tile_time tile_anomaly

%figure;
figure(1); clf

yyaxis right
plot(tile_time,tile_anomaly,'-','linewidth',2)
ylim([-2.5 2.5]);
ylabel('Fit Residual in K');

yyaxis left
plot(mei_time,mei2,'+-');
xlim([datetime(2002,1,1) datetime(2022,1,1)]);
ylabel('MEI Index');
grid;

xlim([datetime(2002,9,1) datetime(2022,9,1)])

hl = legend('MEI Index','1231 cm^{-1} Residual','Location','best');

