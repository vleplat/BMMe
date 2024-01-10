function evaluation = evalf(mu,B1,A,H,lambda1,L_phi1,Wtilde,k)
%evalf evaluates the function f(\mu_k)=\sum_j W_{jk} - 1

% Size computation
[m,~] = size(B1);
[~,n] = size(H);
% vector of ones
e_m1 = ones(m,1);
e_n1 = ones(n,1);
% Computation of the function
B2_vec = e_m1*(e_n1'*H(k,:)') + lambda1*(A(:,k) - L_phi1*Wtilde(:,k) + e_m1*mu');
W_vec = 1/2*(-B2_vec + (B2_vec.^2 + 4*lambda1*L_phi1*B1(:,k)).^(1/2));

evaluation = sum(W_vec) - 1;
end