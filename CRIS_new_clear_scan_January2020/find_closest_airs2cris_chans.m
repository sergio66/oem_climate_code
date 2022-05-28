function chx = find_closest_airs2cris_chans(ch)

%% see closest_airs2645_to_cris1305.m
load closest_airs2645_to_cris1305_chan

chx1 = closest_cris1305_to_airs2645_chan(ch);

%plot(delta_closest_cris1305_to_airs2645_chan,'.-')
closeby = find(abs(delta_closest_cris1305_to_airs2645_chan) < 0.5);
chx = intersect(closeby,chx1);

