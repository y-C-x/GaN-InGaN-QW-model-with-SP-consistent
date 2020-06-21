function psiNorm = findWaveFx(H0,E,n)
% n = 1 for HH
% n = 2 for LH
H = H0-E*eye(length(H0));

psi = null(H);

psiNorm = psi/sqrt(sum(psi.^2));
psiNorm = psiNorm(n:3:end);
end

