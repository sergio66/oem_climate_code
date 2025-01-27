#grep -inR plot_ecmwf_or_era_16days_tile_timestep ../MATLABCODE     >& ugh1find
#sed '/No such file or directory/d' ugh1find

#grep -inR plot_ecmwf_or_era_16days_tile_timestep ../MATLABCODE_Git >& ugh2find
#sed '/No such file or directory/d' ugh2find

find /home/sergio/MATLABCODE/oem_pkg_run/ -type f -name "*.m" -print0 | xargs -I {} -0  grep -in "all64_anomflux" "{}"  \; print  >& ugh3find
sed '/No such file or directory/d' ugh3find

watch "ls -lt ugh*find"
