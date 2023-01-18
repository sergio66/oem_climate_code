iWgtVers = 3;

if iWgtVers == 1
  wgtA = 1; wgtB = 0;
  if driver.iLat >= iLatX & driver.iLat <= 64-iLatX     
     wgtA = 1; wgtB = 0;
  else
     wgtA = 0; wgtB = 1;
  end

elseif iWgtVers == 2
  wgtA = 1; wgtB = 0;
  if driver.iLat >= iLatX+1 & driver.iLat <= 64-(iLatX+1)
     fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,2);
     wgtA = 1; wgtB = 0;
  elseif driver.iLat <= iLatX-1 | driver.iLat >= 64-(iLatX-1)
     wgtA = 0; wgtB = 1;
     fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,0);
  else
     if driver.iLat <= 32
       wgtA = (driver.iLat+1  - iLatX)/2; wgtB = 1-wgtA;
       fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,-1);
     elseif driver.iLat > 32
       wgtB = (driver.iLat+1 - (64-iLatX))/2; wgtA = 1-wgtB;
       fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,+1);
     end
  end

elseif iWgtVers == 3
  wgtA = 1; wgtB = 0;
  if driver.iLat >= iLatX+2 & driver.iLat <= 64-(iLatX+2)
     fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,2);
     wgtA = 1; wgtB = 0;
  elseif driver.iLat <= iLatX-2 | driver.iLat >= 64-(iLatX-2)
     wgtA = 0; wgtB = 1;
     fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,0);
  else
     if driver.iLat <= 32
       wgtA = (driver.iLat+1  - iLatX)/3; wgtB = 1-wgtA;
       fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,-1);
     elseif driver.iLat > 32
       wgtB = (driver.iLat+1 - (64-iLatX))/3; wgtA = 1-wgtB;
       fprintf(1,'driver.iLat = %2i %2i \n',driver.iLat,+1);
     end
  end

end

%{
iLatX = 11;
for ii = 1 : 64
  driver.iLat = ii;
  find_wgtA_wgtB;
  mooA(ii) = wgtA; mooB(ii) = wgtB;
end
%}
