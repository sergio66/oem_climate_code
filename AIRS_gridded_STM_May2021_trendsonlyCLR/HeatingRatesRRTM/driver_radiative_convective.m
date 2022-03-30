function [tau,Tout] = driver_radiative_convective(N,Pin,Tin,TSurfIn)

%% see https://brian-rose.github.io/ClimateLaboratoryBook/courseware/sympy-greenhouse.html

Q = 341.3;        %% solar insoltion at TOA W/m2
Qrefl = 101.9;    %% reflected at TOA by clouds and atmophere
alpha = Qrefl/Q;
sigma = 5.67e-8;  %% Stefan Boltzmann

U = zeros(N+1,N+1);
D = zeros(N+1,N+1);

tau = 0.2;

U(1,1) = 1;
for ii = 2 : N+1
  for jj = 1 : ii
    if ii == jj
      U(ii,jj) = tau;
    elseif jj == 1
      U(ii,jj) = (1-tau)^(ii-1);
    else
      U(ii,jj) = (1-tau)^(ii-jj-1);
    end
  end
end

D(N+1,N+1) = 0;
for ii = N : -1 : 1
  for jj = ii : N
    if ii == jj
      D(ii,jj) = tau;
    elseif jj == 1
      D(ii,jj) = (1-tau)^(ii-1);
    else
      D(ii,jj) = (1-tau)^(ii-jj-1);
    end
  end
end
  
if N <= 6
  U
  D
end



