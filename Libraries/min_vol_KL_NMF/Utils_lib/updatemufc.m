function [mu]=updatemufc(Phi,Eta,mu,epsi,rho)

[K,T]=size(Phi);
JN1  = ones(T,1);
doLoop = true;
while doLoop
    mu_prev=mu;
    Mat=(((Phi).^2+8*(mu*JN1').*Eta).^(1/2)-Phi)./(4*mu*JN1'+eps);
    xi=(sum(Mat.^2,2)-rho*ones(K,1));
    Matp=2*Mat.*((2*(8*(mu*JN1').*Eta)./(((Phi).^2+8*(mu*JN1').*Eta).^(1/2))-4*(((Phi).^2+8*(mu*JN1').*Eta).^(1/2)-Phi))./(4*mu*JN1'+eps).^2);
    xip=sum(Matp,2);
    mu=mu-xi./(xip);
    if(max(abs(mu-mu_prev))<=epsi)
        doLoop=false;
    end
end

%flag=1; %uncomment for debugging and check the convergence rate of mu
% figure;
% semilogy(max(xi_save,0)')

end%EOF