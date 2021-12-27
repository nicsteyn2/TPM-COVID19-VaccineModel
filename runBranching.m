function [D, tEnd, elim] = runBranching(P, D, v, tMax, nMax)

% runBranching.m - branching process implementation
%
% Inputs:
%    P - Parameter structure obtained from getPar.m
%    D - Initial case table describing seed cases (see example below)
%    v - 16x1 array where v(i) gives proportion of age-group i vaccinated
%    tMax - maximum time in days before simulation terminates
%    nMax - maximum number of cases before simulation terminates
%
% Outputs:
%    D - Resulting case table
%    tEnd - time in days when simulation terminated
%    elim - is equal to 1 if simulation ended due to no active cases
%
% Example initial case table input:
%    D = table;
%    D.Exp = 0;
%    D.AgeGroup = 5;
%    D.IsVax = 0;
%
% Author: Nic Steyn, Michael Plank
% Affiliation: Te Pūnaha Matatini
% Email: nicholas.steyn@auckland.ac.nz
% Last revision: 26-12-2021


% ------------- Construct clinical next generation matrix --------------

NGMclin = P.u * P.tI * diag(P.ui) * P.C';


% ------------- Append additional case data --------------

D.ID = (1:height(D))'; % Specify intial case IDs
D.Parent = NaN(height(D),1); % Initial cases do not have parents

D.Clin = rand(height(D),1) < P.IDR(D.AgeGroup); % Determine whether cases are clinical
D.Hosp = rand(height(D),1) < (1 - D.IsVax * P.VEd) .* (P.IHR(D.AgeGroup)./P.IDR(D.AgeGroup)) .* D.Clin; % Determine whether hospitalised
D.Died = rand(height(D),1) < (P.IFR(D.AgeGroup)./P.IHR(D.AgeGroup)) .* D.Hosp; % Determine whether died or not
D.tOns = D.Exp + gamrnd(P.tOnsA, P.tOnsB, height(D), 1); D.tOns(D.Clin==0) = inf; % Set symptom onset time

Ri = sum(NGMclin(:, D.AgeGroup))' .* (1 - (1-P.tau)*(1-D.Clin)) .* (1 - D.IsVax * P.VEt); % Define mean of individual reproduction numbers
D.Ri = gamrnd(P.ssk, Ri/P.ssk); % Assign gamma-distributed individual reproduction numbers

D.Tested = rand(height(D),1) < P.pDetectPre .* D.Clin; % Determine if case is detected by symptomatic testing
D.tTested = D.tOns + exprnd(P.tDetect, height(D), 1); D.tTested(D.Tested==0) = inf; % Set detection time if case is tested
D.Traced = zeros(height(D),1); % Initial cases cannot be contact traced
D.tTraced = inf(height(D),1); % So they have a tTraced time of infinited
D.Detect = D.Tested | D.Traced; % Case is detected if either tested or traced
D.tDetect = min(D.tTested, D.tTraced); % And detection time is the earlier of the two


%------------- Run Branching Process --------------

t = 1; % Start on day 1 (and generate infections that occured in time 0 -> 1)
Dactive = D(t < D.Exp + P.maxInfectTime,:); % Get all cases considered active

% While there are active cases, total infections are less than nMax, and time is less than tMax, keep going
while ( (height(Dactive)>0) && (height(D) <= nMax) && (t <= tMax) )
   
    % On the earliest day after detection, contact trace and test all other cases
    if (t == ceil(min(D.tDetect)))
        notDet = (D.tDetect > min(D.tDetect)); nNotDet = sum(notDet); % Get all cases that have not yet been detected
        D.Traced(notDet) = rand(nNotDet,1) < P.pTrace; % Determine if detected by contact tracing
        D.Tested(notDet) = rand(nNotDet,1) < P.pDetectPost .* D.Clin(notDet); % Determine if detected by additional symptomatic testing
        D.tTraced(notDet) = max(D.tOns(notDet), min(D.tDetect)) + exprnd(P.tTrace, nNotDet, 1); D.tTraced(D.Traced==0) = inf; % Set time of contact tracing
        D.tTested(notDet) = max(D.Exp(notDet), min(D.tDetect)) + exprnd(P.tDetect, nNotDet, 1); D.tTested(D.Tested==0) = inf; % Set time of testing
        D.Detect(notDet) = D.Traced(notDet) | D.Tested(notDet);
        D.tDetect(notDet) = min(D.tTraced(notDet), D.tTested(notDet));
    end
    
    if t > min(D.tDetect)
        Lam = P.postDetectionReduc * (t <= Dactive.tDetect) .* P.auc(t - Dactive.Exp) .* Dactive.Ri; % Set individual infection forces post-detection
    else
        Lam = P.auc(t - Dactive.Exp) .* Dactive.Ri; % Set individual infection forces pre-detection
    end
    nChild = poissrnd(Lam); % Get number of new children per active case
    nChildTotal = sum(nChild); % Get total number of new children
    
    if nChildTotal > 0
        
        Dchild = table; % Create table for child cases
        Dchild.Exp = t*ones(nChildTotal, 1); % All child cases exposed on current day
        Dchild.AgeGroup = zeros(sum(nChild),1); % Pre-allocate age-groups
        parentAges = repelem(Dactive.AgeGroup, nChild, 1); % Find parent age-groups
        for ii = 1:nChildTotal % For each child case, randomly draw their age-group according to the parents NGM column
            Dchild.AgeGroup(ii) = randsample(16, 1, true, NGMclin(:,parentAges(ii)));
        end
        Dchild.IsVax = rand(nChildTotal, 1) < v(Dchild.AgeGroup); % Determine if child cases are vaccinated
        Dchild.ID = ((max(D.ID)+1):(max(D.ID)+nChildTotal))'; % Assign child cases IDs
        Dchild.Parent = repelem(Dactive.ID, nChild, 1); % Assign child cases their parent
        Dchild.Clin = rand(nChildTotal, 1) < P.IDR(Dchild.AgeGroup); % Determine if child cases are clinical
        
        Ri = (1 - Dchild.IsVax * P.VEt) .* sum(NGMclin(:, Dchild.AgeGroup))' .* (1 - (1-P.tau)*(1-Dchild.Clin)); % Determine individual mean reproduction numbers
        Dchild.Ri = gamrnd(P.ssk,  Ri/P.ssk); % Assign gamma distributed individual reproduction numbers
        
        Dchild.Hosp = rand(nChildTotal, 1) < (1 - Dchild.IsVax * P.VEd).*(P.IHR(Dchild.AgeGroup)./P.IDR(Dchild.AgeGroup)) .* Dchild.Clin; % Determine if cases are hospitalised
        Dchild.Died = rand(nChildTotal, 1) < (P.IFR(Dchild.AgeGroup)./P.IHR(Dchild.AgeGroup)) .* Dchild.Hosp; % and if they die
        Dchild.tOns = Dchild.Exp + gamrnd(P.tOnsA, P.tOnsB, nChildTotal, 1); Dchild.tOns(Dchild.Clin==0) = inf; % and assign onset dates to those that are clinical
        
        if t > min(D.tDetect) % Run contact tracing if t > time of first case detection
            Dchild.Tested = rand(nChildTotal, 1) < (P.pDetectPost .* Dchild.Clin); % Determine if case detected by symptomatic testing
            Dchild.tTested = Dchild.tOns + exprnd(P.tDetect, nChildTotal, 1); Dchild.tTested(Dchild.Tested==0) = inf; % Set detection time if tested
            Dchild.Traced = rand(nChildTotal, 1) < P.pTrace; % Determine if case is contact traced
            Dchild.tTraced = Dchild.Exp + exprnd(P.tTrace, nChildTotal, 1); Dchild.tTraced(Dchild.Traced==0) = inf; % Set trace time if traced
            Dchild.Detect = Dchild.Tested | Dchild.Traced; % Case is detected if either tested or traced
            Dchild.tDetect = min(Dchild.tTested, Dchild.tTraced); % And detection time is the earlier of the two
        else % Else don't perform contact tracing
            Dchild.Tested = rand(nChildTotal, 1) < (P.pDetectPre .* Dchild.Clin); % Determine if case detected by symptomatic testing
            Dchild.tTested = Dchild.tOns + exprnd(P.tDetect, nChildTotal, 1); Dchild.tTested(Dchild.Tested==0) = inf; % Set detection time if tested
            Dchild.Traced = zeros(nChildTotal,1); % No contact tracing occurs before detection
            Dchild.tTraced = inf(nChildTotal,1); % So these cases (currently) have a tTraced time of infinity
            Dchild.Detect = Dchild.Tested | Dchild.Traced; % Case is detected if either tested or traced
            Dchild.tDetect = min(Dchild.tTested, Dchild.tTraced); % And detection time is the earlier of the two
        end
        
        keep = rand(nChildTotal, 1) < (1 - P.VEi) * Dchild.IsVax; % Thin child cases according to vax status
        keep(Dchild.IsVax==0) = 1;
        
        D = [D; Dchild(keep==1,:)]; % Append new cases to main case table

    end
    
    t = t + 1; % Increment time
    Dactive = D(t < (D.Exp + P.maxInfectTime),:); % Find new active cases
end

elim = (height(Dactive) == 0); % If the model ended as there were no active cases then it eliminated before tMax/nMax
tEnd = t - 1; % Because time is incremented at the end of the loop, the last day is t - 1


% Append the number of children each case had, and a Reff estimate accounting for censoring (comment out to speed up a little when running multiple iterations)
D.nChild = NaN(height(D),1);
D.Reff = NaN(height(D), 1);
for ii = 1:height(D)
    D.nChild(ii) = sum(D.Parent == D.ID(ii));
    D.Reff(ii) = D.nChild(ii)/wblcdf(tEnd - D.Exp(ii), P.genA, P.genB);
end



end





