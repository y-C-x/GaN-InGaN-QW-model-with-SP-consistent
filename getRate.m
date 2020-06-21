%% Transition Energy
ES = ES.getTransEnergy();
%% OverLap
waveOverLap = overLap(ES);

%% Momentum Matrix Element
MMElement = MME(Fermi,ES,InGaN,waveOverLap); % MB^2 [kg-eV]

%% Gain Spectrum with broadening
sp_rate = emissionRate(waveOverLap,MMElement,wavLen);
sp_rate.nr = nr;
sp_rate = sp_rate.updateC();
sp_rate = sp_rate.getRho(Fermi,ES);
sp_rate = sp_rate.getFermi(ES,Fermi);
sp_rate = sp_rate.getGain();

%% Gain Spectrum with broadening
sp_rate = sp_rate.getGain_broad(gamma,ES);

%% Spontaneous Emission Rate
sp_rate = sp_rate.getSponRate();