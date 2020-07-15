%% Information about different objects
% energyState: Store information about energy states, wave functions and 
% transiton energy between different states.
%
% fermi_level: injection carrier density, effective masses, Fermi levels, 
% carrier distribution
%
% GaN/InN/InGaN_str: material parameters
%
% structure: structure information, band edge profile (Evz, Ecz, Evz0, 
% Ecz0), boundary index, well length, barrier length, energy offset
%
% FDM_C/_H/_P: FDM matrix for conduction band, valence band and Poisson
% equation
%
% strain_eff/sp_pe: strain effect parameters and piezo potential profile
%
% MME: momentum matrix element
%
% overLap: wave overlap profile
%
% sponEmisson: spontaneous emission rate, gain, gain with broadening
%% Code Version: 2020.7.15 - Chengxin
clear
clc
warning('off')
addpath('library')
addpath('functions')
addpath('PreRunData')

%% Constant and material parameters
constant % initialize some global constant
GaN = GaN_str(); % this object contains all parameters for GaN
InN = InN_str(); % this object contains all parameters for InN
x = 0.28; % define the In percentage
Cp = 1.4; % define the bowling parameter
InGaN = InGaN_str(GaN,InN,x,Cp); % this object contains all parameters for InGaN

%% Strain effect
strain = strain_eff(InGaN,GaN); % this object contains info about strain effect

%% Structure Definition
Lb = 103e-10; % barrier length
Lw = 33e-10; % well length
VBO = 0.3; % valence band offset
str = structure(Lb,Lw); % this object contains all info about the structure
str = str.setEoff(VBO,GaN,InGaN); % add the valence band offset
str = str.initialEcv(GaN,InGaN); % initialize the Ec and Ev energy bands
str = str.addStrain(GaN,InGaN,strain); % add strain effect to the structure

%% Sp and Piezo
Pz = sp_pe(GaN,InGaN,str); % this object contains all info about the Pz and Sp polarization
str = str.addPz(Pz); % add the polariztion potential to the structure
str = str.saveOrig(); % save a copy of the current structure profile

%% FDM
% Process of solving the FDM matrix is in FDM.m file. Make change of kt and
% injection carrier density for different cases. I've run this part of code
% already for you so you can just run the code once at 1E17 injection
% carrier density and then load the .mat data provided for the rest part of
% the code.

kt = 0;
n = 1e17 * 1e6; % cm^-3 to m^-3; % injection carrier density
threshold = 0;
% FDM;

%% OR Load Data from PreRunData Directly
% these data is generated under 10.3/3.3 nm GaN/InGaN QW with In of 28%

load('at1E18.mat')
load('at1E19.mat')
% ...
% you can load data you want included in the folder PreRunData or specified
% the injection carrier density you want in the previous section and run
% the FDM agian by uncommenting the FDM command.


%% Get Rate
% I move part of the code to the file getRate.m. It calculates the
% transition energy and wave overlap, as well as emission rate.

wavLen = 400:0.01:650; % nm
nr = 2.7756; % refractive index
gamma = 30e-3; % broadening
getRate;

%% Plotting
hold on
% plot(sp_rate.wavLen*1e9,sp_rate.g_linewidth) % gain  [1/m]
% plot(sp_rate.wavLen*1e9,sp_rate.g_linewidth_broad) % gain with broadening [1/m]
% plot((sp_rate.wavLen*1e9),(sp_rate.r_sp)) % spontaneous emission rate [1/eV-s-m^3]
% plot(sp_rate.hw,sp_rate.r_sp/1e6)
plot(sp_rate.wavLen*1e9,sp_rate.r_sp_broad/1e6)
set(gca, 'XDir','reverse')
% plot(sp_rate.hw,sp_rate.g_linewidth_broad/100) % gain with broadening [1/cm]
% title(sprintf('Total spontaneous emission rate [n = %.1e cm^{-3}]',n/1e6));
grid on
xlabel('wave length [nm]')
% xlabel('hw [eV]')
% xlim([2.7 3.2])
% ylabel('gain [1/cm]');
ylabel('[1/eV-s-cm^3]')
% xlim([550 600]);
% ylim([-1e3 1e4])