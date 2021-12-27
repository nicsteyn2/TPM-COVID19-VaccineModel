function P = getPar(r0)
% getPar - Loads data and sets default parameter values for vaccine models
%
% Inputs:
%   r0 - Reproduction number before any vaccine/infection induced immunity
%
% Outputs:
%    P - parameter structure P
%
% Structure can be edited later using dot-indexing. If changing the
% R0 value (e.g. P.R0 = x) make sure to update the NGM too.
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 27-12-2021


%------------- Fundamental Parameters --------------

P.R0 = r0; % R0
P.tau = 0.5; % Relative infectiousness of subclinicals
P.VEi = 0.7; % Default vaccine efficacy against infection
P.VEt = 0.5; % Default vaccine efficacy against transmission given infection
P.VEd = 0.8; % Default vaccine efficacy against severe disease (and death) given infection

%------------- SEIR Parameters --------------

P.tE = 2.55; % Time spent in exposed compartment
P.tI = 5; % Time spent in infectious compartment
P.tH = 8; % Time spent in hospital

%------------- Branching Process Parameters --------------

P.ssk = 0.5; % Overdispersion/superspreading parameter k
P.maxInfectTime = 21; % Maximum infectious perod (for computational efficiency)
P.genA = 5.67; P.genB = 2.83; % Generation time distribution parameters
P.tOnsA = 5.8; P.tOnsB = 0.95; % Exposure to Onset distribution parameters

P.tDetect = 4; % Mean time from onset to detection & isolation of symptomatic cases
P.pDetectPre = 0.12; P.pDetectPost = 0.4; % Probability of detecting symptomatic case
P.tTrace = 6; % Mean time from exposure to isolation from contact tracing
P.pTrace = 0.7; % Probability of detecting a case by contact tracing

P.postDetectionReduc = 1; % Relative reduction in Reff. 1 implies no reduction.
P.auc = wblcdf((1:P.maxInfectTime)', P.genA, P.genB) - wblcdf((1:P.maxInfectTime)' - 1, P.genA, P.genB);

%------------- Load and Fit Disease Rate Data --------------

nzPopDist = readmatrix("data/nzpopdist.xlsx"); % Load NZ population structure from data folder
fraserIDR = 1 - [0.456; 0.412; 0.370; 0.332; 0.296; 0.265; 0.238; 0.214; 0.192]; % Fraser group symptomatic rates
verityIHR = [0; 0.000408; 0.0104; 0.0343; 0.0425; 0.0816; 0.118; 0.166; 0.184]; % Verity et al hospitalisation rates
verityIFR = [0.0000161; 0.0000695; 0.000309; 0.000844; 0.00164; 0.00595; 0.0193; 0.0428; 0.078]; % Verity et al fatality rates
daviesSusc = [0.4 0.38 0.79 0.86 0.8 0.82 0.88 0.74]; % Davies et al age-structured infectiousness
dataAges = [5; 15; 25; 35; 45; 55; 65; 75; 85]; % Age groups of data

finalAges = (2.5:5:77.5)'; % Define final age groups (matching contact matrix)
P.IDR = interp1(dataAges, fraserIDR, finalAges); % Interpolate within-group estimates
P.IDR(1) = fraserIDR(1); P.IDR(16) = (P.IDR(16)*nzPopDist(16,2) + fraserIDR(9)*sum(nzPopDist(17:end,2)))/sum(nzPopDist(16:end,2));

P.IHR = interp1(dataAges, verityIHR, finalAges); % Interpolate within-group estimates
P.IHR(1) = verityIHR(1); P.IHR(16) = (P.IHR(16)*nzPopDist(16,2) + verityIHR(9)*sum(nzPopDist(17:end, 2)))/sum(nzPopDist(16:end, 2)); % Fit out-of-group estimates

P.IFR = interp1(dataAges, verityIFR, finalAges); % Interpolate within-group estimates 
P.IFR(1) = verityIFR(1); P.IFR(16) = (P.IFR(16)*nzPopDist(16,2) + verityIFR(9)*sum(nzPopDist(17:end, 2)))/sum(nzPopDist(16:end, 2)); % Fit out-of-group estimates

P.ui = interp1(dataAges(1:8), daviesSusc, finalAges); % Interpolate within-group estimates
P.ui(1) = daviesSusc(1); P.ui(end) = daviesSusc(end);

P.IHR(1) = 1e-26; % Set to near zero to avoid division by zero error in ODE
P.IFR(1) = 0;

%------------- Specify Population Structure --------------

P.totalPopSize = 5e6; % Total size of assumed population
popDist = zeros(16,1); % Create popDist vector
popDist(1:15) = nzPopDist(1:15,2); % Fill entries with population distribution
popDist(16) = sum(nzPopDist(16:end,2)); % Aggregate 75+ age-groups
P.popDist = popDist/sum(popDist); % Save to structure
P.age = finalAges; % Specify corresponding age-groups
P.popCount = P.popDist * P.totalPopSize; % Specify popDist in terms of implied group sizes

%------------- Load Contact Matrix and Define NGM --------------

C = readmatrix("data/nzcontmatrix.xlsx"); % Get Prem et al contact matrix from data folder
P.C = zeros(16,16); % Create our contact matrix
for ii = 1:16
    for jj = 1:16
        P.C(ii,jj) = 0.5*(C(ii,jj) + (P.popDist(jj)/P.popDist(ii)) * C(jj,ii)); % Force detailed balance condition
    end
end
NGM = P.tI*diag(P.ui)*P.C'*diag(P.IDR + P.tau*(1-P.IDR)); % Construct "first guess at the NGM"
P.u = P.R0/max(abs(eig(NGM))); % Choose u to give desired R0
P.NGM = P.u*NGM; % Set final NGM







