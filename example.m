% example.m - some examples that are easier to implement than the results
% in the paper
%
% Use the run section button to run code one block at a time
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 26-12-2021



%% Example 1: SEIR -------------------------------------------------------
% Other m-files required: getPar.m, runSEIR.m

clc
clear
close all

% Inputs:
r0 = 4.5; % Set the basic reproduction number (default 4.5)
propVax = zeros(16,1); propVax(14:16) = 0.8; % Set proportion of each age group vaccinated (default 80% of over 65yrs)
dailyImports = 5; % Set daily imported cases (default to 5 per-day)
tMax = 365; % Set maximum time in days (default to 1 year)

% Setup:
P = getPar(r0); %Load parameter structure (see getPar.m)

ic = zeros(16,16); % Create initial conditions matrix
ic(:,1) = (1 - propVax) .* P.popDist; % Assign individuals to the unvaccinated class
ic(:,6) = propVax .* P.popDist; % And the rest to the vaccinated class

m = dailyImports*P.popDist/sum(P.popDist); % Set daily imported cases (default spread proportionally over all ages)

% Run model:
D = runSEIR(P, ic, m, tMax); % Run SEIR model

% Aggregate SEIR curves
S = zeros(tMax,1); EIA = zeros(tMax, 1); R = zeros(tMax, 1); F = zeros(tMax, 1);

for t = 1:tMax
    S(t) = sum(D.S(t,:) + D.Sv(t,:));
    EIA(t) = sum(D.E(t,:) + D.Ev(t,:) + D.I(t,:) + D.Iv(t,:) + D.A(t,:) + D.Av(t,:));
    F(t) = sum(D.F(t,:) + D.Fv(t,:));
    R(t) = sum(D.R(t,:) + D.Rv(t,:)) - F(t);
end

plot(1:t, [S EIA R F], 'LineWidth', 2)
xlim([0 100])
legend(["Susceptible", "Exposed/Infected", "Recovered", "Died"], 'Location', "East")
xlabel("Time (days)"); ylabel("Proportion of population")




%% Example 2: Branching Process -------------------------------------------
% The function call getTracingReduc(P) on line 77 uses a parfor loop. If
% this is unavailable then leave P.tTrace, P.pTrace, P.tDetect, and
% P.pDetectPost at their default values. In this case getTracingReduc(P)
% returns the value 0.4370.
%
% Other m-files required: getPar.m, runBranching.m, support/getTracingReduc.m


clc
clear
close all
addpath('support')


% Inputs
r0 = 4.5; % Set the basic reproduction number (default 4.5)
propVax = zeros(16,1); propVax(14:16) = 0.8; % Set proportion of each age group vaccinated (default 80% of over 65yrs)
tMax = inf; nMax = inf; % Set no stopping criteria (will only stop if eliminated)
postDetectReduc = 0.2; % Set basic reproduction number, before tracing/isolation, to x% of r0 after outbreak detection 


% Setup:
P = getPar(r0); % Load parameter structure (see getPar.m)
P.postDetectionReduc = postDetectReduc; % Assign postDetectionReduc in parameter structure

D0 = table; D0.Exp = 0; D0.AgeGroup = 5; D0.IsVax = 0; % Initialise by infecting an unvaccinated 25-29 y/o


% Check post-detection R < 1 (else model will never stop, see line 65)
fprintf("Calculating post-detection R...\n")
postDetectR = P.postDetectionReduc * getTracingReduc(P) * r0;
if postDetectR > 0.99
    error("Elimination not guaranteed, model will run forever if continued.")
else
    fprintf("Rt = %.2f (after first detection)\nAll set. Beginning model...\n", postDetectR)
end

%% Run model
% This code is separated so can run multiple times without re-calculating
% post-detection R.

[D, tEnd, elim] = runBranching(P, D0, propVax, tMax, nMax);

clc
fprintf("Total cases: %d\n", height(D))
if sum(D.Detect) == 0
    fprintf("Outbreak eliminated prior to detection\n")
else
    fprintf("Time to first detection: %d days\n", floor(min(D.tDetect)))
    fprintf("Time from first detection until elimination: %d days\n", tEnd - floor(min(D.tDetect)))
    fprintf("Total hospitalisations: %d\n", sum(D.Hosp))
    fprintf("Total fatalities: %d\n", sum(D.Died))
end




%% Example 3: Reff Calculation --------------------------------------------
% Other m-files required: support/getReff.m, support/getVaccineRollout.m

clc
clear
close all
addpath("support")

% Inputs:
r0 = 3.0; % Basic reproduction number
nVax = (0:1e4:4.5e6)'; % Number of vaccine doses to allocate
maxVax = 0.9; % Vaccinate up to maxVax of each age-group

% Setup:
P = getPar(r0); % Load default parameter structure
[V, ageChanges] = getVaccineRollout(P, nVax, maxVax); % Get default vaccine rollout with doses given by nVax and up to maxVax of each age-group vaccinated

% Calculate Reff:
Reff = getReff(P,V);

% Find when Reff < 1
ind = find(Reff < 1, 1);

if isempty(ind)
    fprintf("Either Reff < 1 always, or vaccination with maxVax = %.2f is insufficient to get Reff < 1.\n", maxVax)
else
    fprintf("Reff < 1 once %d vaccine courses have been administered (%.0f%% of the popn).\n", nVax(ind), sum(V(:,ind) .* P.popDist)*100);
    fprintf("%.0f%% of 65+ y/o, %.0f%% of 15-64 y/o, %.0f%% of 0-14 y/o.\n", V(16,ind)*100, V(13, ind)*100, V(3,ind)*100)
end





