clear variables;
close all;
clc;



t=linspace(0,5,100);
sig=sin(t);
sig=sig';
xn=zeros(1,100);
gamma=0.001;
lambda=[0.01,0.1,0.5,1,10,50];
x_final=zeros(6,100);

fonction_cout=zeros(6,1);



H = matH(size(sig),'gaussian',7);

z = H*sig+0.5*randn(length(sig),1);

operateur='laplacien';


switch operateur
    
    case 'identité'
        
        I = matGamma(size(sig),'identity');
        
        for i=1:length(lambda)

            xn1=xn'-2*gamma*lambda(i)*(I')*I*(xn')-gamma*2*(H')*(H*(xn')-(z));
            eps=10^(-6);
            iteration=1;


            while(norm(abs(xn1-xn'))>eps)

                xn=xn1';
                xn1=xn'-2*gamma*lambda(i)*(I')*I*(xn')-gamma*2*(H')*(H*(xn')-(z));
                fonction_cout(i,iteration)=norm(H*xn1-z)^2+lambda(i)*norm(I*xn1)^2;
                iteration=iteration+1;
            end
            x_final(i,:)=xn1;

        end
        
    case 'gradient'
        
        G = matGamma(size(sig),'gradient');
        
        for i=1:length(lambda)

            xn1=xn'-2*gamma*lambda(i)*(G')*G*(xn')-gamma*2*(H')*(H*(xn')-(z));
            eps=10^(-6);
            iteration=1;


            while(norm(abs(xn1-xn'))>eps)

                xn=xn1';
                xn1=xn'-2*gamma*lambda(i)*(G')*G*(xn')-gamma*2*(H')*(H*(xn')-(z));
                fonction_cout(i,iteration)=norm(H*xn1-z)^2+lambda(i)*norm(G*xn1)^2;
                iteration=iteration+1;

            end
            x_final(i,:)=xn1;

        end

        
    case 'laplacien'
        
        L = matGamma(size(sig),'laplacian');
        
        for i=1:length(lambda)

            xn1=xn'-gamma*2*(H')*(H*(xn')-(z))-2*gamma*lambda(i)*(L')*L*(xn');
            eps=10^(-6);
            iteration=1;

            while(norm(abs(xn1-xn'))>eps)

                xn=xn1';
                xn1=xn'-gamma*2*(H')*(H*(xn')-(z))-2*gamma*lambda(i)*(L')*L*(xn');
                fonction_cout(i,iteration)=norm(H*xn1-z)^2+lambda(i)*norm(L*xn1)^2;
                iteration=iteration+1;

            end
            x_final(i,:)=xn1;

        end
        
        
end



for i=1:length(lambda)

    figure(1)

    subplot(3,2,i)

    plot(t,x_final(i,:));
    hold on;
    plot(t,z);
    plot(t,sig);

    title(['reconstruction par Tikhonov ',operateur,' lambda = ',num2str(lambda(i)),' et gamma = ',num2str(gamma)]);
    legend('signal recréé','signal observé','signal attendu');
    xlabel('temps');
    ylabel('amplitude');

    figure(2)

    subplot(3,2,i);

    plot(fonction_cout(i,:));
    title(['Fonction de cout ',operateur,' lambda = ',num2str(lambda(i)),' et gamma = ',num2str(gamma)]);
    xlabel('nombre d itération');
    ylabel('valeur de la fonction');
    
end



