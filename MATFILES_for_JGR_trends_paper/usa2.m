function h = usa;

%% usa flag colors
%% see http://www.mikesoltys.com/2011/09/16/tool-of-the-week-prettyer-matlab-plots/

m = 064; dark = 10; light= 20;
m = 128; dark = 20; light= 40;


B(1,:)= [zeros(1,dark) linspace(0,1,light) ones(1,light) linspace(1,.5,dark)];
B(2,:)= [zeros(1,dark) linspace(0,1,light) linspace(1,0,light) zeros(1,dark)];
B(3,:)= [linspace(.5,1,dark) ones(1,light) linspace(1,0,light) zeros(1,dark)];

h = (B');
