% this function finds KL divergence D(X,WH)
function dist = KLobj(X,W,H)
   beta=1;
   dist=betaDiv(X+eps,W*H+eps,beta);   

end
