hold on;
run_retrieval_latbins_AIRS_loop
plot(lat,co2,'b+-')
save Test/alldriver_fit_obs_with_robust_error_p001err alldriver
whos alldriver
alldriver
alldriver(33)
alldriver(33).oem
g = load('Test/alldriver_fit_obs');
g.alldriver(33).oem
g
whos
whos g
whos alldriver
whos g
whos Test/alldriver_fit_obs_with_robust_error_p001err alldriver
whos g
1135655942 - 1135655766
g
alldriver(20)
alldriver(20).lls
alldriver(20).jacobian
g.alldriver(20).jacobian
alldriver(20).oem
figure(1)
figure(2)
imagesc(alldriver(20).oem.inv_se)
imagesc(alldriver(20).oem.se)
imagesc(alldriver(20).oem.inv_se)
colorbar
caxis([0 1E4])
caxis([0 1E3])
imagesc(g.alldriver(20).oem.inv_se)
driver.rateset
clf
plot(driver.rateset.unc_rates)
driver
driver.rateset
driver.jacobian
ig = driver.jacobian.chanset;
plot(driver.rateset.unc_rates(ig))
hold on;
plot(g.driver.rateset.unc_rates(ig),'r')
clf
plot(alldriver(18).rateset.unc_rates(ig))
hold on;
plot(g.alldriver(18).rateset.unc_rates(ig),'r')
plot(alldriver(34).rateset.unc_rates(ig),'k')
plot(g.alldriver(34).rateset.unc_rates(ig),'g')
hl = legend('JOB18','JOB18','JOB34','JOB34')
grid
figure(1)
xlim([-80 80])
h1=plot(alldriver(34).rateset.unc_rates(ig),'m')
figure(2)
delete(h1)
h1=plot(alldriver(34).rateset.unc_rates(ig),'m')
delete(h1)
h1=plot(alldriver(35).rateset.unc_rates(ig),'m')
delete(h1)
h1=plot(alldriver(21).rateset.unc_rates(ig),'m')
h2=plot(alldriver(23).rateset.unc_rates(ig),'c')
figure(3)
plot(f(ig),alldriver(23).rateset.unc_rates(ig),'b')
hold on;
plot(f(ig),alldriver(19).rateset.unc_rates(ig),'r')
grid
plot(f(ig),alldriver(32).rateset.unc_rates(ig),'g')
figure(4)
nc_cor = nc_rates(alldriver(19);
nc_cor = nc_rates(alldriver(19));
plot(f,nc_cor)
nc_cor = nc_rates(alldriver(32));
hold on;
plot(f,nc_cor,'r')
nc_cor = nc_rates(alldriver(35));
plot(f,nc_cor,'g')
1.65/1.1
1.45/1.25
clf
nc_cor = nc_rates(alldriver(19));
plot(f,nc_cor.*alldriver(19).rateset.unc_rates,'b-')
plot(f9ig),nc_cor(ig).*alldriver(19).rateset.unc_rates(ig),'b-')
plot(f(ig),nc_cor(ig).*alldriver(19).rateset.unc_rates(ig),'b-')
hold on;
nc_cor = nc_rates(alldriver(32));
plot(f(ig),nc_cor(ig).*alldriver(32).rateset.unc_rates(ig),'r-')
nc_cor = nc_rates(alldriver(35));
plot(f(ig),nc_cor(ig).*alldriver(35).rateset.unc_rates(ig),'g-')
0.014/0.0045
0.016/0.0075
alldriver(1)
alldriver(1).jacobian
k = load(alldriver(1).jacobian.filename);
k
figure(5)
plot(k.f,squeeze(k.M_TS_jac_all(19,:,1)))
hold on;
plot(k.f,squeeze(k.M_TS_jac_all(32,:,1)),'r')
plot(k.f,squeeze(k.M_TS_jac_all(35,:,1)),'g')
grid
figure(6);
find(f > 704,1)
f(190)
plot(squeeze(k.M_TS_jac_all(:,191,1)),'g')
figure(6)
plot(squeeze(k.M_TS_jac_all(:,191,1)),'g')
plot(squeeze(k.M_TS_jac_all(:,191,1)),'b+-')
plot(save_lat,squeeze(k.M_TS_jac_all(:,191,1)),'b+-')
plot(lat,squeeze(k.M_TS_jac_all(:,191,1)),'b+-')
k
plot(lat,squeeze(k.M_TS_jac_all(:,191,170)),'b+-')
plot(lat,squeeze(k.M_TS_jac_all(:,191,1)),'b+-')
hold on;
plot(lat,squeeze(k.M_TS_jac_all(:,191,170)),'r+-')
plot(lat,squeeze(k.M_TS_jac_all(:,191,170))*30,'r+-')
clf
plot(lat,squeeze(k.M_TS_jac_all(:,191,1)),'b+-')
hold on;
plot(lat,-squeeze(k.M_TS_jac_all(:,191,170))*50,'r+-')
grid
12/8
7.5/0.075
0.075/0.075
0.075/0.035
k
cjac = squeeze(k.M_TS_jac_all(:,:,1));
tjac = squeeze(k.M_TS_jac_all(:,:,104:200));
whos cjac
whos tjac
for i=1:36
[u(i,:,:) s(i,:,:) v(i,:,:)]=svd(squeeze(tjac(i,:,:)));
i
end
whos u
whos s
whos v
figure(7)
plot(f,squeeze(u(1,1,:)))
plot(f,squeeze(u(1,end,:)))
plot(f,squeeze(u(1,:,1)))
plot(f,squeeze(u(1,:,2)))
plot(f,squeeze(u(1,:,20)))
plot(f,squeeze(u(1,:,1)))
doc svd
imagesc(squeeze(s(20,:,:)))
clf
imagesc(squeeze(s(20,:,:)))
colorbar
plot(diag(squeeze(s(20,:,:))))
plot(diag(squeeze(s(20,:,:))),'+-')
plot(f,squeeze(u(1,:,1)))
plot(diag(squeeze(s(20,:,:))),'+-')
plot(f,squeeze(u(1,:,1)))
plot(diag(squeeze(s(20,:,:))),'+-')
imagesc(squeeze(s(20,:,:)))
colorbar
caxis([0 0.005])
caxis([0 0.0005])
caxis([0 0.00005])
whos s
s = s(:,1:97,:);
whos s
imagesc(squeeze(s(20,:,:)))
whos s
for i=1:36
diags(i,:) = diag(squeeze(i,:,:));
end
for i=1:36
i
whos diags
diags
for i=1:36
diags(i,:) = diag(squeeze(s(i,:,:)));
end
whos diags
plot(diags(:,1))
plot(diags(:,2))
plot(diags(:,3))
plot(diags(:,4))
plot(diags(:,5))
plot(diags(:,6))
plot(f,squeeze(u(1,:,1)))
dot(squeeze(u(1,:,1)),squeeze(u(1,:,20))
dot(squeeze(u(1,:,1)),squeeze(u(1,:,20)))
dot(squeeze(u(1,:,1)),squeeze(u(1,:,1)))
whos u
plot(f,squeeze(u(1,:,1)))
hold on;
plot(f,cjac,'r')
whos cjac
clf
plot(f,squeeze(u(1,:,1)))
hold on;
plot(f,cjac(1,:),'r')
cc = dot(cjac(1,:),cjac(1,:))
cc = dot(cjac(1,:)/cc,cjac(1,:)/cc)
cc = dot(cjac(1,:)/sqrt(cc),cjac(1,:)/sqrt(cc))
cc = dot(cjac(1,:),cjac(1,:))
type dot
whos cjac
cc = dot(cjac(1,:),cjac(1,:))
dot(cjac(1,:)/(cc*36),cjac(1,:)/(cc*36))
dot(cjac(1,:)/(cc),cjac(1,:)/(cc))
sum(cjac(1,:).^2)
dot(cjac(1,:)/sqrt(cc)),cjac(1,:)/sqrt((cc)))
dot(cjac(1,:)/sqrt(cc)),cjac(1,:)/sqrt((cc))))
dot(cjac(1,:)/sqrt(cc),cjac(1,:)/sqrt(cc))
ncjac1 = cjac(1,:)/sqrt(cc);
dot(squeeze(u(1,:,1)),ncjac1)
dot(squeeze(u(1,:,2)),ncjac1)
dot(squeeze(u(1,:,20)),ncjac1)
dot(squeeze(u(1,:,3)),ncjac1)
dot(squeeze(u(1,:,4)),ncjac1)
dot(squeeze(u(1,:,1)),ncjac1)
for i=1:36
cjac_norm(i,:) = cjac(i,:)/sqrt(dot(cjac(i,:),cjac(i,:));
for i=1:36
cjac_norm(i,:) = cjac(i,:)/sqrt(dot(cjac(i,:),cjac(i,:)));
end
dot(cjac_norm(20,:),cjac_norm(20,:))
for i=1:36
co_lin(i) = dot(cjac_norm(i,:),squeeze(u(i,:,1)));
end
figure(8)
plot(lat,co_lin,'b+-')
figure(1);
hold on;
hx=plot(lat,(10*co_lin),'go-')
delete(hx);hx=plot(lat,(7*co_lin),'go-')
delete(hx);hx=plot(lat,2-(7*co_lin),'go-')
delete(hx);hx=plot(lat,(7*co_lin),'go-')
for i=1:36
co_lin2(i) = dot(cjac_norm(i,:),squeeze(u(i,:,2)));
end
figure(8)
hold on;
plot(lat,co_lin2,'r+-')
whos diags
figure(9)
imagesc(diags)
plot(lat,co_lin2.*diags(:,2),'r+-')
plot(lat,co_lin2.*diags(:,2)','r+-')
grid
for i=1:36
co_lin3(i) = dot(cjac_norm(i,:),squeeze(u(i,:,3)));
end
plot(lat,co_lin3.*diags(:,3)','g+-')
clf
plot(lat,co_lin,'b+-')