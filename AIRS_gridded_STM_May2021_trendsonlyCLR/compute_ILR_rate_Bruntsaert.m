function [deltaILR,ILR] = compute_ILR_rate_Bruntsaert(eair,Tx2m,dWVsurf,dTx2m)

stefan_boltzmann = 5.6700e-08; %% W/m2/K^4
ILR = 1.24*stefan_boltzmann*(eair./Tx2m).^(1/7) .* (Tx2m.^4);

% Eqn 9, see Bk 48 of notes
% Earth Syst. Dynam., 14, 1363–1374, 2023 https://doi.org/10.5194/esd-14-1363-2023
% Understanding variations in downwelling longwave radiation using Brutsaert’s equation
% Yinglin Tian1,2, Deyu Zhong1, Sarosh Alam Ghausi2,3, Guangqian Wang1, and Axel Kleidon2
deltaILR = stefan_boltzmann*(Tx2m .^4) *1.24/7.*(eair./Tx2m).^(1/7).*(dWVsurf - dTx2m./Tx2m);

%% sergio
deltaILR = ILR .* (dWVsurf  + 27 * dTx2m./Tx2m)/7;
