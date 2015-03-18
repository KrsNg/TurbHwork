clear all;
clc;

size = 6;

format long;

K = 0.41;

U_vonK = zeros(size,1);

U_vanD = zeros(size,1);

y = zeros(size,1);


U_vonK(1,1)= 12.24641;
U_vonK(2,1)=12.68913;
U_vonK(3,1)=13.06368;
U_vonK(4,1)=13.3883;
U_vonK(5,1)=13.67476;
U_vonK(6,1)=13.93107;

U_vanD(1,1) = 18.75522;
U_vanD(2,1) = 19.19794;
U_vanD(3,1) = 19.57251;
U_vanD(4,1) = 19.89713;
U_vanD(5,1) = 20.18358;
U_vanD(6,1) = 20.4399;


counter = 0;

% 
for i = 1 : size
   y(i) = 250 + counter;
   counter = counter + 50;
   
   C_vonK(i,1) = U_vonK(i,1) -  (1/K)*log( y(i) );
   
   C_vanD(i,1) = U_vanD(i,1) -  (1/K)*log( y(i) );
   
   
end;
% 
% plot(x,y);

