% resultsBranching.m - reproduce branching process results
%
% Code to reproduce Figure 3
%
% Other m-files required: getTracingReduc.m, getReff.m, getVaccineRollout.m
% Subfunctions: getStochasticPElim(), getOtherResults()
%
% Author: Nic Steyn, Michael Plank
% Affiliation: Te Pūnaha Matatini
% Email: nicholas.steyn@auckland.ac.nz
% Last revision: 27-12-2021


clc
clear
close all
% addpath("support")

% Simulation parameters (full results use nRuns = 1e5)
nRuns = 10;

% Setup
col = [205 4 2; 236 116 0; 155 190 0]/255; % colours for plotting
nVax = (0:0.5e6:4.5e6)'; propVax = nVax/5e6;
V = getVaccineRollout(getPar(3.0), nVax, 0.9);


%% Get all results (this takes a while for high nRuns) --------------------

tic

% R0 = 3.0 ------
P = getPar(3.0);
% Low detection
P.pDetectPre = 0.05;
R3.low.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R3.low.nInf, R3.low.tElim, R3.low.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);
% Moderate detection
P.pDetectPre = 0.1;
R3.med.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R3.med.nInf, R3.med.tElim, R3.med.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);
% High detection
P.pDetectPre = 0.15;
R3.high.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R3.high.nInf, R3.high.tElim, R3.high.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);

% R0 = 4.5 ------
P = getPar(4.5);
% Low detection
P.pDetectPre = 0.05;
R4.low.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R4.low.nInf, R4.low.tElim, R4.low.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);
% Moderate detection
P.pDetectPre = 0.1;
R4.med.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R4.med.nInf, R4.med.tElim, R4.med.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);
% High detection
P.pDetectPre = 0.15;
R4.high.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R4.high.nInf, R4.high.tElim, R4.high.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);

% R0 = 6.0 ------
P = getPar(6.0);
% Low detection
P.pDetectPre = 0.05;
R6.low.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R6.low.nInf, R6.low.tElim, R6.low.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);
% Moderate detection
P.pDetectPre = 0.1;
R6.med.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R6.med.nInf, R6.med.tElim, R6.med.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);
% High detection
P.pDetectPre = 0.15;
R6.high.pElim = getStochasticPElim(P, V, nRuns, 1000);
[R6.high.nInf, R6.high.tElim, R6.high.nHosps] = getOtherResults(P, V, nRuns, inf, 0.25);

x = toc;
fprintf("\n\n\nTime elapsed: %.0fs for nRuns = %d\n", x, nRuns)


%% Plot -------------------------------------------------------------------

% Infections at Detection Plots

ymax = max( [max(max(R3.low.nInf)) max(max(R4.low.nInf)) max(max(R6.low.nInf))] )*1.05;

subplot(4,3,1)
errorbar(propVax-0.01, R3.low.nInf(1,:), R3.low.nInf(1,:)-R3.low.nInf(2,:), R3.low.nInf(3,:)-R3.low.nInf(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R3.med.nInf(1,:), R3.med.nInf(1,:)-R3.med.nInf(2,:), R3.med.nInf(3,:)-R3.med.nInf(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R3.high.nInf(1,:), R3.high.nInf(1,:)-R3.high.nInf(2,:), R3.high.nInf(3,:)-R3.high.nInf(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
ylabel("Infections at Detection");
title(["Basic Reproduction Number: 3.0","(a)"]);
ylim([0 ymax]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,2)
errorbar(propVax-0.01, R4.low.nInf(1,:), R4.low.nInf(1,:)-R4.low.nInf(2,:), R4.low.nInf(3,:)-R4.low.nInf(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R4.med.nInf(1,:), R4.med.nInf(1,:)-R4.med.nInf(2,:), R4.med.nInf(3,:)-R4.med.nInf(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R4.high.nInf(1,:), R4.high.nInf(1,:)-R4.high.nInf(2,:), R4.high.nInf(3,:)-R4.high.nInf(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
title(["Basic Reproduction Number: 4.5","(b)"]);
ylim([0 ymax]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,3)
errorbar(propVax-0.01, R6.low.nInf(1,:), R6.low.nInf(1,:)-R6.low.nInf(2,:), R6.low.nInf(3,:)-R6.low.nInf(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R6.med.nInf(1,:), R6.med.nInf(1,:)-R6.med.nInf(2,:), R6.med.nInf(3,:)-R6.med.nInf(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R6.high.nInf(1,:), R6.high.nInf(1,:)-R6.high.nInf(2,:), R6.high.nInf(3,:)-R6.high.nInf(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
title(["Basic Reproduction Number: 6.0","(c)"]);
ylim([0 ymax]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');


% Probability of elimination

ymin = min( [min(min(R3.high.pElim)), min(min(R4.high.pElim)), min(min(R6.high.pElim))] )*0.95;

subplot(4,3,4)
plot(propVax, R3.low.pElim, '--*', 'Color', col(1,:), 'LineWidth', 1.5, 'MarkerSize', 10);
hold on
plot(propVax, R3.med.pElim, '--*', 'Color', col(2,:), 'LineWidth', 1.5, 'MarkerSize', 10);
plot(propVax, R3.high.pElim, '--*', 'Color', col(3,:), 'LineWidth', 1.5, 'MarkerSize', 10);
ylabel("P(elim) before 1000 cases");
title("(d)");
ylim([ymin 1]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,5)
plot(propVax, R4.low.pElim, '--*', 'Color', col(1,:), 'LineWidth', 1.5, 'MarkerSize', 10);
hold on
plot(propVax, R4.med.pElim, '--*', 'Color', col(2,:), 'LineWidth', 1.5, 'MarkerSize', 10);
plot(propVax, R4.high.pElim, '--*', 'Color', col(3,:), 'LineWidth', 1.5, 'MarkerSize', 10);
title("(e)");
ylim([ymin 1]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,6)
plot(propVax, R6.low.pElim, '--*', 'Color', col(1,:), 'LineWidth', 1.5, 'MarkerSize', 10);
hold on
plot(propVax, R6.med.pElim, '--*', 'Color', col(2,:), 'LineWidth', 1.5, 'MarkerSize', 10);
plot(propVax, R6.high.pElim, '--*', 'Color', col(3,:), 'LineWidth', 1.5, 'MarkerSize', 10);
title("(f)");
ylim([ymin 1]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');


% Time to elimination

ymax = max( [max(max(R3.low.tElim)), max(max(R4.low.tElim)), max(max(R6.low.tElim))] )*1.05;

subplot(4,3,7)
errorbar(propVax-0.01, R3.low.tElim(1,:), R3.low.tElim(1,:)-R3.low.tElim(2,:), R3.low.tElim(3,:)-R3.low.tElim(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R3.med.tElim(1,:), R3.med.tElim(1,:)-R3.med.tElim(2,:), R3.med.tElim(3,:)-R3.med.tElim(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R3.high.tElim(1,:), R3.high.tElim(1,:)-R3.high.tElim(2,:), R3.high.tElim(3,:)-R3.high.tElim(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
ylabel(["Time to Elimination", "(days)"]);
title("(g)");
ylim([0 ymax]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,8)
errorbar(propVax-0.01, R4.low.tElim(1,:), R4.low.tElim(1,:)-R4.low.tElim(2,:), R4.low.tElim(3,:)-R4.low.tElim(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R4.med.tElim(1,:), R4.med.tElim(1,:)-R4.med.tElim(2,:), R4.med.tElim(3,:)-R4.med.tElim(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R4.high.tElim(1,:), R4.high.tElim(1,:)-R4.high.tElim(2,:), R4.high.tElim(3,:)-R4.high.tElim(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
title("(h)");
ylim([0 ymax]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,9)
errorbar(propVax-0.01, R6.low.tElim(1,:), R6.low.tElim(1,:)-R6.low.tElim(2,:), R6.low.tElim(3,:)-R6.low.tElim(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R6.med.tElim(1,:), R6.med.tElim(1,:)-R6.med.tElim(2,:), R6.med.tElim(3,:)-R6.med.tElim(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R6.high.tElim(1,:), R6.high.tElim(1,:)-R6.high.tElim(2,:), R6.high.tElim(3,:)-R6.high.tElim(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
title("(i)");
ylim([0 ymax]); xlim([-0.01 0.91]); xticks(propVax);
xline(0.7296,'HandleVisibility','off');


% Total hospitalisations

ymax = max( [max(max(R3.low.nHosps)), max(max(R4.low.nHosps)), max(max(R6.low.tElim))] )*1.05;

subplot(4,3,10)
errorbar(propVax-0.01, R3.low.nHosps(1,:), R3.low.nHosps(1,:)-R3.low.nHosps(2,:), R3.low.nHosps(3,:)-R3.low.nHosps(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R3.med.nHosps(1,:), R3.med.nHosps(1,:)-R3.med.nHosps(2,:), R3.med.nHosps(3,:)-R3.med.nHosps(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R3.high.nHosps(1,:), R3.high.nHosps(1,:)-R3.high.nHosps(2,:), R3.high.nHosps(3,:)-R3.high.nHosps(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
ylabel("Total Hospitalisations");
title("(j)");
xlim([-0.01 0.91]); xticks(propVax);
% ylim([0 ymax]);
ylim([0 max(max(R3.low.nHosps))*1.05]);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,11)
errorbar(propVax-0.01, R4.low.nHosps(1,:), R4.low.nHosps(1,:)-R4.low.nHosps(2,:), R4.low.nHosps(3,:)-R4.low.nHosps(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R4.med.nHosps(1,:), R4.med.nHosps(1,:)-R4.med.nHosps(2,:), R4.med.nHosps(3,:)-R4.med.nHosps(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R4.high.nHosps(1,:), R4.high.nHosps(1,:)-R4.high.nHosps(2,:), R4.high.nHosps(3,:)-R4.high.nHosps(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
xlabel("Proportion Fully Vaccinated"); 
title("(k)");
xlim([-0.01 0.91]); xticks(propVax);
% ylim([0 ymax]); 
ylim([0 max(max(R4.low.nHosps))*1.05]);
xline(0.7296,'HandleVisibility','off');

subplot(4,3,12)
errorbar(propVax-0.01, R6.low.nHosps(1,:), R6.low.nHosps(1,:)-R6.low.nHosps(2,:), R6.low.nHosps(3,:)-R6.low.nHosps(1,:), '-x', 'Color', col(1,:), 'LineWidth', 3);
hold on
errorbar(propVax, R6.med.nHosps(1,:), R6.med.nHosps(1,:)-R6.med.nHosps(2,:), R6.med.nHosps(3,:)-R6.med.nHosps(1,:), '-x', 'Color', col(2,:), 'LineWidth', 3);
errorbar(propVax+0.01, R6.high.nHosps(1,:), R6.high.nHosps(1,:)-R6.high.nHosps(2,:), R6.high.nHosps(3,:)-R6.high.nHosps(1,:), '-x', 'Color', col(3,:), 'LineWidth', 3);
title("(l)");
xlim([-0.01 0.91]); xticks(propVax);
ylim([0 max(max(R6.low.nHosps))*1.05]); % ylim([0 ymax]); 
xline(0.7296,'HandleVisibility','off');




%% Sub-functions ----------------------------------------------------------

function pElim = getStochasticPElim(P, V, nRuns, nMax)
% Estimates the probability of elimination

fprintf("Running (R0 = %.1f, pDetectPre = %.2f): pElim.\n", P.R0, P.pDetectPre);

szV = size(V);

didElim = zeros(nRuns, szV(2));

D0 = table;
D0.Exp = 0;
D0.IsVax = 0;
D0.AgeGroup = 0;

P.postDetectionReduc = 1;

for jj = 1:szV(2) % For each vaccine allocation
    v = V(:,jj);
    parfor ii = 1:nRuns % run the model nRuns times
        D = D0; D.AgeGroup = randsample(16,1,true,P.popDist);
        [~, ~, elim] = runBranching(P, D, v, inf, nMax);
        didElim(ii,jj) = elim; % and count how many eliminated
    end
end

pElim = mean(didElim);

end


function [nInf, tElim, nHosps] = getOtherResults(P, V, nRuns, nMax, alpha)
% Get other results from stochastic sims

fprintf("Running (R0 = %.1f, pDetectPre = %.2f): other.\n", P.R0, P.pDetectPre);

szV = size(V);
infAtDet = zeros(nRuns, szV(2)); tToElim = infAtDet; totalHosps = infAtDet;
D0 = table; D0.Exp = 0; D0.IsVax = 0; D0.AgeGroup = 0;

traceReduc = getTracingReduc(P);

for jj = 1:szV(2)
    
    P.postDetectionReduc = 1/3;
    v = V(:,jj);
    
    if ((1 - traceReduc) * P.postDetectionReduc * getReff(P,V(:,jj)) < 0.98)% If R with control is < 0.98 and there is testing run the full results
        parfor ii = 1:nRuns
            % Randomly choose age-group and run model
            D = D0; D.AgeGroup = randsample(16,1,true,P.popDist);
            [D, tEnd, elim] = runBranching(P, D, v, inf, nMax);
            
            % Save results
            infAtDet(ii,jj) = sum(D.Exp < min(D.tDetect));
            if sum(D.Detect(D.Exp <= tEnd))==0; infAtDet(ii,jj) = NaN; end
            
            tToElim(ii,jj) = max(D.Exp) - min(D.tDetect) + 14;
            if (elim == 0) || sum(D.Detect)==0; tToElim(ii,jj) = NaN; end
            
            totalHosps(ii,jj) = sum(D.Hosp);
            if sum(D.Detect)==0; totalHosps(ii,jj) = NaN; end
        end
    else % Otherwise only get size at detection and leave others blank. Setting postdetectreduc to 0 forces the sim to stop (else would run forever)
        P.postDetectionReduc = 0;
        parfor ii = 1:nRuns
            % Randomly choose age-group and run model
            D = D0; D.AgeGroup = randsample(16,1,true,P.popDist);
            [D, tEnd, ~] = runBranching(P, D, v, inf, nMax);
            
            % Save results
            infAtDet(ii,jj) = sum(D.Exp < min(D.tDetect));
            if sum(D.Detect(D.Exp <= tEnd))==0; infAtDet(ii,jj) = NaN; end 

            tToElim(ii,jj) = NaN;
            totalHosps(ii,jj) = NaN;
        end
    end
        
    
end

nInf = [nanmedian(infAtDet); quantile(infAtDet, alpha); quantile(infAtDet, 1-alpha)];
tElim = [nanmedian(tToElim); quantile(tToElim, alpha); quantile(tToElim, 1-alpha)];
nHosps = [nanmedian(totalHosps); quantile(totalHosps, alpha); quantile(totalHosps, 1-alpha)];

end

