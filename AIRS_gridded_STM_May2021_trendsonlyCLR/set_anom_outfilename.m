if driver.i16daytimestep > 0 & topts.ocb_set == 2 & driver.NorD > 0
  if iInd <= 5000
    zanom_outdir = ['OutputAnomaly/00001_05000/'];
  elseif iInd <= 10000
    zanom_outdir = ['OutputAnomaly/05001_10000/'];
  elseif iInd <= 15000
    zanom_outdir = ['OutputAnomaly/10001_15000/'];
  elseif iInd <= 20000
    zanom_outdir = ['OutputAnomaly/15001_20000/'];
  elseif iInd <= 25000
    zanom_outdir = ['OutputAnomaly/20001_25000/'];
  elseif iInd <= 30000
    zanom_outdir = ['OutputAnomaly/25001_30000/'];
  elseif iInd <= 35000
    zanom_outdir = ['OutputAnomaly/30001_35000/'];
  elseif iInd <= 40000
    zanom_outdir = ['OutputAnomaly/35001_40000/'];
  end
elseif driver.i16daytimestep > 0 & topts.ocb_set == 2 & driver.NorD < 0
  if iInd <= 5000
    zanom_outdir = ['OutputAnomaly_DAY/00001_05000/'];
  elseif iInd <= 10000
    zanom_outdir = ['OutputAnomaly_DAY/05001_10000/'];
  elseif iInd <= 15000
    zanom_outdir = ['OutputAnomaly_DAY/10001_15000/'];
  elseif iInd <= 20000
    zanom_outdir = ['OutputAnomaly_DAY/15001_20000/'];
  elseif iInd <= 25000
    zanom_outdir = ['OutputAnomaly_DAY/20001_25000/'];
  elseif iInd <= 30000
    zanom_outdir = ['OutputAnomaly_DAY/25001_30000/'];
  elseif iInd <= 35000
    zanom_outdir = ['OutputAnomaly_DAY/30001_35000/'];
  elseif iInd <= 40000
    zanom_outdir = ['OutputAnomaly_DAY/35001_40000/'];
  end
elseif driver.i16daytimestep > 0 & (topts.ocb_set == 0 | topts.ocb_set == 1) & driver.NorD > 0
  zanom_outdir = ['OutputAnomaly_CAL/35001_40000/'];    
  error('kjgkljlsjg')
elseif driver.i16daytimestep > 0 & (topts.ocb_set == 0 | topts.ocb_set == 1) & driver.NorD < 0
  zanom_outdir = ['OutputAnomaly_CAL_DAY/35001_40000/'];    
  error('kjgkljlsjg')
end
