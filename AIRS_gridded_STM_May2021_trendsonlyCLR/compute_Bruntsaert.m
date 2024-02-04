function [ILR,ewater,emiss_atmosphere] = compute_Bruntsaert(tdew2m,tair2m)

% Earth Syst. Dynam., 14, 1363–1374, 2023 https://doi.org/10.5194/esd-14-1363-2023
% Understanding variations in downwelling longwave radiation using Brutsaert’s equation
% Yinglin Tian1,2, Deyu Zhong1, Sarosh Alam Ghausi2,3, Guangqian Wang1, and Axel Kleidon2

ewater = 6.1079 *exp(17.269.*(tdew2m-273)./(237.3+(tdew2m-273)));  %% in mb, Eqn 7 of paper

emiss_atmosphere = 1.24 * (ewater./tair2m).^(1/7);                 %% no units

stefan_boltzmann = 5.6700e-08; %% W/m2/K^4

%ILR = 1.24*stefan_boltzmann*(ewater./tair2m).^(1/7) .* (tair2m.^4);  %% W/m2
ILR = stefan_boltzmann * emiss_atmosphere .* (tair2m.^4);             %% W/m2
