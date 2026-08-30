filename = '../../oem_pkg_run_sergio_AuxJacs/Aux_jacs_Generic/SARTA_CLOUD_O3_JACS/sarta_M_TS_jac_all_5_4_87_97_97_cld.mat';
varname  = 'M_TS_jac_all';
scalar_i = 1:5;
scalar_i = 1:9;
water_i  = 10:106;
temp_i   = 107:203;
ozone_i  = 204:300;

load(filename);

for ii = 1 : 5
  figure(ii); clf; colormap jet
end

iLat = input('enter latbin (1-40) : ');
while iLat > 0 & iLat < 41
  figure(1); clf
    pcolor(f,1:97,squeeze(M_TS_jac_all(iLat,:,temp_i))'); shading flat; colorbar; set(gca,'ydir','reverse'); title('T(z) jac')
  figure(2); clf
    pcolor(f,1:97,squeeze(M_TS_jac_all(iLat,:,water_i))'); shading flat; colorbar; set(gca,'ydir','reverse'); title('WV(z) jac')
  figure(3); clf
    pcolor(f,1:97,squeeze(M_TS_jac_all(iLat,:,ozone_i))'); shading flat; colorbar; set(gca,'ydir','reverse'); title('O3(z) jac')
  figure(4);
    plot(f,squeeze(M_TS_jac_all(iLat,:,1:5))); title('col jacs and stemp')
  figure(5);
    plot(f,squeeze(M_TS_jac_all(iLat,:,6:9))); title('cld jac')
  iLat = input('enter latbin (1-40) : ');    
end
