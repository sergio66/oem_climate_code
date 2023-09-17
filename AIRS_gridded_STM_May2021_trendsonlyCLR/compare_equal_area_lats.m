boo = equal_area_spherical_bands(32);
plot(0.5*(boo(1:end-1)+boo(2:end)),diff(boo),'.-'); xlim([-90 +90]);

whos rlat
plot(0.5*(boo(1:end-1)+boo(2:end)),diff(boo),'.-',rlat,diff(rlat65)); xlim([-90 +90]);
plot(0.5*(boo(1:end-1)+boo(2:end)),diff(boo),'.-',rlat,diff(rlat65),'.'); xlim([-90 +90]);
plot(0.5*(boo(1:end-1)+boo(2:end)),diff(boo),'.-',rlat,diff(rlat65),'.-'); xlim([-90 +90]);
