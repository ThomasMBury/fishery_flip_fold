

% MATLAB script to produce simulated data for the exploitation and agetruncation
% scenarios under environmental and process error
% Estimate CV, AR1 and Nonlinearity DeltaRho and produce basic plots
% 28 April 2016 by V Dakos
clear all
close all
% Stochasticity strength
sr0 = 0.05; % demeographic noise
se0 = 0.1; % environmental noise
su0 = 0; % observational noise
% Simulation settings
iter0 = 1000; % number of realisations
steps0 = 1000; % number of increments in control parameter
T0 = 1000; % number of time0-steps
% age truncation scenario
rstart = 0.5;
rend = 2.6; % stop before choas : hopf at r=2
r_values = linspace(rstart,rend,steps0); % vector of control parameter
% exploitation scenario
Fstart = 0;
Fend = 3;
F_values = linspace(Fstart,Fend,steps0); % vector of control parameter
% simulated data
X_agetruncation = Ricker_agetruncation(sr0,se0,su0,iter0,r_values,T0);
X_exploitation = Ricker_exploitation(sr0,se0,su0,iter0,F_values,T0);


