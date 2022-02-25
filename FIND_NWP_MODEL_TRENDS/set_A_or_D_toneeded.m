%Airs_Temp    = (Airs_Temp_A + Airs_Temp_D)/2;
%Airs_STemp   = (Airs_STemp_A + Airs_STemp_D)/2;
%Airs_H2OVap  = (Airs_H2OVap_A + Airs_H2OVap_D)/2;
%Airs_Ozone   = (Airs_Ozone_A + Airs_Ozone_D)/2;
%Airs_OLR     = (Airs_OLR_A + Airs_OLR_D)/2;
%Airs_Clr_OLR = (Airs_ClrOLR_A + Airs_ClrOLR_D)/2;

if iDorA > 0
  Airs_Temp    = Airs_Temp_D;
  Airs_STemp   = Airs_STemp_D;
  Airs_H2OVap  = Airs_H2OVap_D;
  if iL3orCLIMCAPS == +1
    Airs_RHSurf  = Airs_RHSurf_D;
  else
    Airs_RHSurf = 0 * Airs_STemp;
  end
  Airs_RH      = Airs_RH_D;
  Airs_Ozone   = Airs_Ozone_D;
  Airs_CO      = Airs_CO_D;
  Airs_CH4     = Airs_CH4_D;
  Airs_CldPres  = Airs_CldPres_D;
  if iL3orCLIMCAPS == +1
    Airs_OLR     = Airs_OLR_D;
    Airs_ClrOLR = Airs_ClrOLR_D;
    
    Airs_LiqWater = Airs_LiqWater_D;
    Airs_IceT     = Airs_IceT_D;
    Airs_IceSze   = Airs_IceSze_D;
    Airs_IceOD    = Airs_IceOD_D;
    Airs_CldFrac  = Airs_CldFrac_D;
  else
    Airs_OLR = ones(size(Airs_STemp));
    Airs_ClrOLR = ones(size(Airs_OLR));
    Airs_LiqWater = ones(size(Airs_OLR));
    Airs_IceT = ones(size(Airs_OLR));
    Airs_IceSze = ones(size(Airs_OLR));
    Airs_IceOD = ones(size(Airs_OLR));
    Airs_CldFrac = ones(size(Airs_OLR));
    Airs_CldFrac = Airs_Temp_D(:,1,:,:); Airs_CldFrac = ones(size(Airs_CldFrac));
  end
else
  Airs_Temp    = Airs_Temp_A;
  Airs_STemp   = Airs_STemp_A;
  Airs_H2OVap  = Airs_H2OVap_A;
  if iL3orCLIMCAPS == +1
    Airs_RHSurf  = Airs_RHSurf_A;
  else
    Airs_RHSurf = 0 * Airs_STemp;
  end
  Airs_RH      = Airs_RH_A;
  Airs_Ozone   = Airs_Ozone_A;
  Airs_CO      = Airs_CO_A;
  Airs_CH4     = Airs_CH4_A;
  Airs_CldPres  = Airs_CldPres_A;

  if iL3orCLIMCAPS == +1
    Airs_OLR     = Airs_OLR_A;
    Airs_ClrOLR = Airs_ClrOLR_A;
    
    Airs_LiqWater = Airs_LiqWater_A;
    Airs_IceT     = Airs_IceT_A;
    Airs_IceSze   = Airs_IceSze_A;
    Airs_IceOD    = Airs_IceOD_A;
    Airs_CldFrac  = Airs_CldFrac_A;
  else
    Airs_OLR = ones(size(Airs_STemp));
    Airs_ClrOLR = ones(size(Airs_OLR));
    Airs_LiqWater = ones(size(Airs_OLR));
    Airs_IceT = ones(size(Airs_OLR));
    Airs_IceSze = ones(size(Airs_OLR));
    Airs_IceOD = ones(size(Airs_OLR));
    Airs_CldFrac = Airs_Temp_A(:,1,:,:); Airs_CldFrac = ones(size(Airs_CldFrac));
  end
end
