clear variables;
close all;
clc;


t=linspace(0,5,100);
z=sin(t);
xn=zeros(1,100);
gamma=0.05;

plot(t,xn);
hold on;
plot(t,z);

H = matH(size(z),'gaussian',7);

xn1=xn'-gamma*2*(H')*(H*(xn')-(z'));
eps=10^(-1);
iteration=0;

while(norm(abs(xn1-xn'))>eps)
    
    xn=xn1';
    xn1=xn'-gamma*2*(H')*(H*(xn')-(z'));
    iteration=iteration+1;

end
plot(t,xn1)

disp(iteration);



