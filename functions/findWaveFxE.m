function psiNorm = findWaveFxE(H0,E)
H = H0-E*eye(length(H0));

psi = null(H);

psiNorm = psi/sqrt(sum(psi.^2));
end
