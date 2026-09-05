%  oni = load('ONI_sep2023.txt');
%  [aaa,bbb] = size(oni);
%  oniS = oni(1,1); oniE = oni(aaa,1);
%  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
%  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;

% oniX = 'ONI_sep2023.txt';
oni_data = 'ONI_sep2023.txt';
oni_data = 'ONI_jul2026.txt';  
if exist(oni_data)
  read_oni
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
  plot(onidd,oni); title('ONI Oceanic Nino Index')
end
