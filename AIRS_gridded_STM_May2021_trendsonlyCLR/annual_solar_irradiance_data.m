function solar = annual_solar_irradiance_data()

%% https://www.ncei.noaa.gov/data/total-solar-irradiance/access/yearly/
%% wget https://www.ncei.noaa.gov/data/total-solar-irradiance/access/yearly/tsi_v02r01_yearly_s1610_e2023_c20240221.nc

%{
ncdump says time(time) ;
       	    time:units = "days since 1610-01-01 00:00:00

	    float TSI(time) ;
	    	  TSI:long_name = "NOAA Climate Data Record of Yearly Total Solar Irradiance (W m-2)" ;
                  TSI:standard_name = "solar_irradiance" ;
		  TSI:units = "W m-2" ;
%}

solar = read_netcdf_lls('tsi_v02r01_yearly_s1610_e2023_c20240221.nc');
solar.solartime = 1610.0 + solar.time/365.25;
solar.solartime = 1610.5 + ((1:414)-1);
plot(solar.solartime,solar.TSI); xlim([1700 2025]); grid
plot(solar.solartime,solar.TSI); xlim([2002 2025]); grid
