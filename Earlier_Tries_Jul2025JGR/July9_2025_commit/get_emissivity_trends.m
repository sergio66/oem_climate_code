function emiseffect = get_emissivity_trends(driver)

%emistrends = load('/home/sergio/MATLABCODE/CAMEL_emissivity/Trends_camelV003_emis4608tiles/emistrendsV003.mat');
emistrends = load('/home/sergio/MATLABCODE/CAMEL_emissivity/Trends_camelV003_emis4608tiles/emistrendsV003.mat','t00','t0X');
emiseffect = emistrends.t0X(:,driver.iibin) - emistrends.t00(:,driver.iibin);

