%% FDM
while (threshold == 0) % normally the program should reach convergence within 20 loops
m_h = FDM_H(str,GaN); % this object contains all info for solving the k.p 6x6 Hamiltonian 
m_h = m_h.initialFDM(str,strain,GaN,InGaN); 
m_h = m_h.buildMatrix(kt); % construct the matrix at kt = 0;

m_c = FDM_C(str,GaN); % this object contains all info for solving the Schroodinger Eq
m_c = m_c.initialFDM(str,strain,GaN,InGaN);
m_c = m_c.buildMatrix(kt); % construct the matrix at kt = 0;

%% Energy State and Wave Function
ES = energyState(m_h,m_c); % this object contains all info about energy states and wave functions
ES = ES.getWave(m_h,m_c); % calculate the wave function for every energy states

%% Quasi Fermi Level & Charge Distribution
Fermi = fermi_level(InGaN,strain,n); % this object contains all info about the fermi level and DOS
Fermi = Fermi.solveFermi(str,ES); % solve for the fermi level

%% Poisson
Poisson = FDM_P(str,GaN); % this object contains all info for solving the poisson equation
Poisson = Poisson.initialFDM(str,InGaN);
Poisson = Poisson.buildMatrix();
Poisson = Poisson.solveVsc(Fermi.rho); % solve for the added potential

%% Update Str
str = str.updateStr(Poisson.Vsc); % add the potential to the original structure
%% Convergence test
fprintf('%.4e eV\n',max(abs(Poisson.Vsc)));
if max(abs(Poisson.Vsc)) < 0.001
    fprintf('Current Precision: %.4e eV\n',max(abs(Poisson.Vsc)))
    threshold = 1;
    break;
end
end