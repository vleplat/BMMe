function c = bisectionMethod(a,b,epsilon,B1,A,H,lambda1,L_phi1,Wtilde,k) %(ensure change of sign between a and b) error=1e-4

% Sanity check
if (evalf(a,B1,A,H,lambda1,L_phi1,Wtilde,k)*evalf(b,B1,A,H,lambda1,L_phi1,Wtilde,k))>=0
    disp('no change of signs in f')
    return
end

maxIter = 10^4;
% Bisection method
c=a+(b-a)/2;
i = 0;
while abs(evalf(c,B1,A,H,lambda1,L_phi1,Wtilde,k))>epsilon && abs(b-a)>epsilon && i<=maxIter
    if (evalf(c,B1,A,H,lambda1,L_phi1,Wtilde,k)*evalf(a,B1,A,H,lambda1,L_phi1,Wtilde,k))<0
        b=c;
    else
        a=c;
    end
    c=a+(b-a)/2;
    i = i + 1;
end

if (i==maxIter)
    warning("maximum iteration reached for Bisection method, check the results !")
end
