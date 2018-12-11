
function X = Ricker_agetruncation(sr0,se0,su0,iter0,r_values,T0)
% RICKER_EXPLOITATION function to simulate timeseries to estimate AR, VAR and DeltaRho
% 28 April 2016 by V Dakos
% Initialise random number generator
rng('shuffle','twister')
% Simulation settings
T = 1000; % timesteps that are equal to the size of timeseries used
iter = iter0; % number of stochastic realisations
steps = length(r_values); % steps for changing control parameter
initT = 100; % transient time (discarded - simulation starts at equilibirum)
% Parameters
r = 0.75; % growth rate
K = 10; % carrying capacity
F = 0; % grazing rate
P = 2; % exponent for sigmoid functional response
H = 0.75; % half-saturation rate
% Stochasticity strength
sr = sr0; % process error on r
se = se0; % environmental noise
su = su0; % measurement error
% Solving Difference Equation
X = zeros(T,iter);
counter = 1;



extra = 0;
k = 1;
while k <= iter
x = K;
for t = 1:initT+T
    r = r_values(1);
    if t>initT
        r=r_values(ceil((t-initT)*(steps/T))); % increase r incrementally after transient
    end % increase F incrementally after transient
b = r/K;
r1 = r*exp(sr*randn(1)); % process error
x = (x*exp(r1-b*x)-F*x^P/(H^P+x^P))*exp(se*randn(1));
if x<0
extra = extra+1;
break
end
if t>initT % store values after discarding transient
X(t-initT,counter) = x*exp(su*randn(1)); % observation error
end
end
if x<0 % if negative solutions, discard and rerun realisation
k = k;
counter = counter;
else
k=k+1;
counter = counter+1;
end

end
filename = sprintf('sim_data/Ricker_changer_sr%g_se%g.csv',sr,se);
csvwrite(filename,X)


