%{
ON CHEYENNE
THIS SQUISHES DOWN THE *h0* monthly files from 1.9 G to 87 M by saving of only the required variables
## make sure you do    module load nco/5.0.7
for file in /glade/scratch/aherring/archive/cam.clubbmf.bsort_FHIST_ne30pg2_ne30pg2_mg17_L58dev_3600pes_220720_run76_Nx2yrs/atm/hist/*.cam.h0.*; do
   echo $file
   ncks -v lat,lon,lev,time,date,co2vmr,ch4vmr,n2ovmr,CLDTOT,CLOUD,CLDICE,CLDLIQ,IWC,LANDFRAC,TS,T,Q,O3,PS,U10 $file -o /glade/scratch/sergio/adam_1998_2022/reduced_$(basename $file)
done
%}

%{
see https://project.cgd.ucar.edu/projects/CLUBB-MF/runs/ for description of runs 91/92
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE
fin0 = '/umbc/xfs2/strow/asl/s1/sergio/JUNK/cam.clubbmf.bsort_FHIST_ne30pg2_ne30pg2_mg17_L58dev_3600pes_220720_run76_Nx2yrs.cam.h0.2010-03.nc';
fin = '/umbc/xfs2/strow/asl/s1/sergio/JUNK/reduced_cam.clubbmf.bsort_FHIST_ne30pg2_ne30pg2_mg17_L58dev_3600pes_220720_run76_Nx2yrs.cam.h0.2010-03.nc';
fout = '/umbc/xfs2/strow/asl/s1/sergio/JUNK//cam_h0_2010-03.mat';

c0 = read_netcdf_lls(fin0);
c = read_netcdf_lls(fin);

p.nlevs = 58;
p.lat = c.lat;
p.lon = c.lon;
p.plevs = c.lev;
p.time = c.time;
p.date = c.date;
p.co2vmr = c.co2vmr;
p.ch4vmr = c.ch4vmr;
p.n2ovmr = c.n2ovmr;

p.tcc = c.CLDTOT;
p.cc   = c.CLOUD;
p.ciwc = c.CLDICE;
p.clwc = c.CLDLIQ;
p.iwc  = c.IWC;
p.lwc  = c.LWC;

p.landfrac = c.LANDFRAC;

p.stemp = c.TS;
p.ptemp = c.T;
p.gas_1 = c.Q;
p.gas_3 = c.O3;
p.spres = c.PS;

p.u10 = c.U10;
p.v10 = c.V10;

saver = ['save ' fout ' p'];
if ~exist(fout)
  eval(saver)
  eval(['!ls -lth ' fin '  ' fout])
end
