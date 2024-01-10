function [mu]=updatemu(C,S,D,W,mu,epsi)

[F,K]=size(W);
JF1=ones(F,1);
doLoop = true;
while doLoop
    mu_prev=mu;
    Mat=W .*(((C+JF1*mu').^2+S).^(1/2)-(C+JF1*mu'))./(D+eps);
    xi=(sum(Mat,1)-ones(1,K))';
    Matp=(W./D+eps).*((C+JF1*mu')./(((C+JF1*mu').^2+S).^(1/2))-ones(F,K));
    xip=sum(Matp,1)';
    mu=mu-xi./xip;
    if(max(abs(mu-mu_prev))<=epsi)
        doLoop=false;
    end
end

% flag=1; %uncomment for debugging and check the convergence rate of mu
% figure;
% semilogy(max(xi_save,0)')

end%EOF