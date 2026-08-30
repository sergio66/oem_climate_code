[ppmvLAY1,ppmvAVG1,ppmvMAX1,pavgLAY1,tavgLAY1,ppmv500_1,ppmv75_1,ppmvSURF_1] = layers2ppmv(h,p,1:length(p.stemp),1);

duO3_0 = dobson_rtp(h,p);
duO3_F = dobson_rtp(h,pert);

duO3_col0   = dobson_gas_rtp(h, p, 3);
duO3_300mb0 = dobson_gas_rtp(h, p, 3, 300);
[ppmvLAY3,ppmvAVG3,ppmvMAX3,pavgLAY3,tavgLAY3,ppmv500_3,ppmv75_3,ppmvSURF_3] = layers2ppmv(h,p,1:length(p.stemp),3);

duO3_colpert   = dobson_gas_rtp(h, pert, 3);
duO3_300mbpert = dobson_gas_rtp(h, pert, 3, 300);
[ppmvLAYpert3,ppmvAVGpert3,ppmvMAXpert3,pavgLAYpert3,tavgLAYpert3,ppmv500pert3,ppmv75pert3,ppmvSURFpert3] = layers2ppmv(h,pert,1:length(p.stemp),3);
[ppmvLAYpert_unc3,ppmvAVGpert_unc3,ppmvMAXpert_unc3,pavgLAYpert_unc3,tavgLAYpert_unc3,ppmv500pert_unc3,ppmv75pert_unc3,ppmvSURFpert_unc3] = layers2ppmv(h,pert_unc,1:length(p.stemp),3);

[nlayO3,~] = size(ppmvLAY3);

fracO3 = pert.gas_3 ./ p.gas_3 - 1;
fracO3 = fracO3 .* (ones(101,1) * maskLF);

fracO3unc = pert_unc.gas_3 ./ p.gas_3 - 1;
fracO3unc = pert.gas_3_unc ./ p.gas_3;

deltaO3 = ppmvLAYpert3 - ppmvLAY3;
deltaO3 = deltaO3 .* (ones(nlayO3,1) * maskLF);

deltaO3unc = ppmvLAYpert_unc3 - ppmvLAY3;
deltaO3unc = deltaO3unc .* (ones(nlayO3,1) * maskLF);

