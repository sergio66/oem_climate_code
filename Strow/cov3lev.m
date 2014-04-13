function y = cov3lev(c);

x1 = 1:30;
x2 = 31:97;

pt1 = 1:c.trans1+ceil(c.trans2-c.trans1)/2;
pt2 = pt1(end)+1:97;

m1 = (c.lev2-c.lev1)/2;
m2 = (c.lev3-c.lev2)/2;

y1 = m1*tanh(c.width1*(pt1-c.trans1)) + (c.lev2 + c.lev1)/2;
y2 = m2*tanh(c.width2*(pt2-c.trans2)) + (c.lev3 + c.lev2)/2;

keyboard

y = [y1 y2];

% Sample values for c structure
%
% c.trans1 = 20;
% c.trans2 = 73;
% c.width1 = 1/2;
% c.width2 = 1/3;
% c.lev1 = 1;
% c.lev2 = 3;
% c.lev3 = 4;

