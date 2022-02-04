clear variables;
close all;
clc;



xn=18;
gamma=0.3;

xn1=xn-gamma*2*xn;
eps=10^(-6);
iteration=0;

while(abs(xn1-xn)>eps)
    xn=xn1;
    xn1=xn-gamma*2*xn;
    iteration=iteration+1;
end

disp(xn1);
disp(iteration);