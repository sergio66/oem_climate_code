if ia_OorC_DataSet_Quantile(1) == 2 
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
    fprintf(1,'input_spectrum_number : locallatbin localtimestep ---> procnumber = %4i : %4i %4i -----> %4i \n',input_spectrum_number,latnumber,iLocalTimeStep,procnumber);
    mapAnomData_to_processor(input_spectrum_number,1) = latnumber;
    mapAnomData_to_processor(input_spectrum_number,2) = iLocalTimeStep;
    mapAnomData_to_processor(input_spectrum_number,3) = procnumber;
  end
  clear latnumber iLocalTimeStep procnumber iPreviouslatnumber input_spectrum_number

  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
end
