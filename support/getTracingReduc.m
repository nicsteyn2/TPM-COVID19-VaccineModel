function reduc = getTracingReduc(P)

% getTracingReduc.m - calculate reduction in reproduction number due to
% contract tracing & case isolation. See supplementary section 1.2 "Case
% detection, isolation and controls".
%
% Inputs:
%    P - Parameter structure obtained from getPar.m
%
% Outputs:
%    reduc - reduction in R due to contact tracing and isolation
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 26-12-2021

N = 1000000; % Size of stochastic simulation

rhoTrace = zeros(N,1); % P(Tg > Ttrace)
rhoDetect = zeros(N,1); % P(Tg > Tons + Tdetect)
rhoBoth = zeros(N,1); % P(Tg > min(Ttrace, Tons + Tdetect))

parfor ii = 1:N
    Tg = ceil(wblrnd(P.genA, P.genB)); % Sample Tg
    tTrace = exprnd(P.tTrace); % Sample Ttrace
    tDetect = exprnd(P.tDetect) + gamrnd(P.tOnsA, P.tOnsB); % Sample Tons + Tdetect
    
    rhoTrace(ii) = tTrace < Tg;
    rhoDetect(ii) = tDetect  < Tg;
    rhoBoth(ii) = min(tTrace, tDetect) < Tg;
end

rhoTrace = mean(rhoTrace);
rhoDetect = mean(rhoDetect);
rhoBoth = mean(rhoBoth);

% Calculate updated NGM
NGM = P.NGM;
for jj = 1:16
    NGM(:,jj) = (1 - P.pTrace*(1-P.pDetectPost*P.IDR(jj))*rhoTrace - (1-P.pTrace)*P.pDetectPost*P.IDR(jj)*rhoDetect - P.pTrace*P.pDetectPost*P.IDR(jj)*rhoBoth)*NGM(:,jj);
end

% Calculate reduction in Reff compared to baseline
reduc = 1 - max(abs(eig(NGM)))/max(abs(eig(P.NGM)));

end