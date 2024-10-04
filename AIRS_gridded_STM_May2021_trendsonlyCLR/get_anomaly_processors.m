function [jobjunk,mapAnomData_to_processor] = get_anomaly_processors(ia_OorC_DataSet_Quantile,iNumAnomTimeSteps,iNumAnomTiles,iNumAnomJobsPerProc,JOB)

if ia_OorC_DataSet_Quantile(1) == 2 
  junkprevious = -1;
  jobjunk = struct;
  junkcnt = 0;
  disp('I suggest running test_run_retrieval_setlatbin_AIRS_loop_anomaly.m BEFOREHAND to see how many processors you need')
  disp('may be more than your usual 64!')
  clear mapperAnom2Processor
  iPreviouslatnumber = 0;
  for input_spectrum_number = 1 : iNumAnomTimeSteps * iNumAnomTiles
    procnumber = floor((input_spectrum_number-1)/(iNumAnomJobsPerProc) + 1);
    latnumber = floor((input_spectrum_number-1)/(iNumAnomTimeSteps) + 1);
    if latnumber == iPreviouslatnumber
      iLocalTimeStep = iLocalTimeStep  + 1;
    else  
      iLocalTimeStep = 1;
    end
    iPreviouslatnumber = latnumber;
    if procnumber == JOB
      fprintf(1,'input_spectrum_number : locallatbin localtimestep ---> procnumber = %4i : %4i %4i -----> %4i clustjob = %03i \n',input_spectrum_number,latnumber,iLocalTimeStep,procnumber,JOB);
      junkcnt = junkcnt + 1;
      jobjunk.inputspectrum(junkcnt) = input_spectrum_number;
      jobjunk.latnumber(junkcnt) = latnumber;
      jobjunk.localtimestep(junkcnt) = iLocalTimeStep;
    else
      if junkprevious ~= procnumber
        fprintf(1,'input_spectrum_number : locallatbin localtimestep ---> procnumber = %4i : %4i %4i -----> %4i \n',input_spectrum_number,latnumber,iLocalTimeStep,procnumber);
        junkprevious = procnumber;
      end
    end 
    mapAnomData_to_processor(input_spectrum_number,1) = latnumber;
    mapAnomData_to_processor(input_spectrum_number,2) = iLocalTimeStep;
    mapAnomData_to_processor(input_spectrum_number,3) = procnumber;
  end
  clear latnumber iLocalTimeStep procnumber iPreviouslatnumber input_spectrum_number

  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))

  if JOB > 0
    junk = [JOB junkcnt jobjunk.inputspectrum(1) jobjunk.inputspectrum(end) jobjunk.latnumber(1) jobjunk.latnumber(end) jobjunk.localtimestep(1) jobjunk.localtimestep(end)];
    fprintf(1,'for JOB = %3i there are %3i spectra that will be processed, from %5i to %5i, with latnumber %2i to %2i and timestep from %4i to %4i \n',junk);
    fprintf(1,'for JOB = %3i there are %3i spectra that will be processed, from %5i to %5i, with latnumber %2i to %2i and timestep from %4i to %4i \n',junk);
    fprintf(1,'for JOB = %3i there are %3i spectra that will be processed, from %5i to %5i, with latnumber %2i to %2i and timestep from %4i to %4i \n',junk);
  end

end
