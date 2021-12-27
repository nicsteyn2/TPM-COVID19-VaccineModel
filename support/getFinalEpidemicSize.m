function [nInfs, nHosps, nDied, maxHosp, tMaxHosp] = getFinalEpidemicSize(P, V, m, tMax)
% getFinalEpidemicSize.m - runs SEIR model with inputs and extracts final
% epidemic size info
%
% Inputs:
%    P - Parameter structure obtained from getPar.m
%    V - 16xN matrix where V(i,j) gives the proportion of age-group i that
%    is vaccinated at time j
%    m - 16x1 array where m(i) is the daily number of infected arrivals of age i
%    tMax - Maximum time of simulation in days
%
% Outputs:
%    nInfs - a Nx3 matrix where nInfs(i,j) gives the total number of
%    infections when the vaccine rollout is described by V(:,i). j = 1
%    includes all infections, j = 2 includes infections in non-vaccinated
%    individuals, and j = 3 includes infections in vaccinated individuals
%    nHosps, nDied, maxHosp, tMaxHosp - as nInfs
%
% Other m-files required: runSEIR.m
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 27-12-2021

popDist = P.popDist;
totalPopSize = P.totalPopSize;
VEi = P.VEi;

% Pre-allocate arrays
szV = size(V);
nInfs = zeros(szV(2), 3);
nHosps = zeros(szV(2), 3);
nDied = zeros(szV(2), 3);
maxHosp = zeros(szV(2), 3);
tMaxHosp = zeros(szV(2), 3);

for ii = 1:szV(2)
    ic = zeros(16,16);
    ic(:,1) = (1 - V(:,ii)) .* popDist;
    ic(:,6) = (1 - VEi) * V(:,ii) .* popDist;
    
    [D, t] = runSEIR(P, ic, m, tMax);
    
    % Save total results
    nInfs(ii,1) = sum(D.R(end,:)+D.Rv(end,:))*totalPopSize;
    nHosps(ii,1) = sum(D.F(end,:)+D.Fv(end,:)+D.Dis(end,:)+D.Disv(end,:))*totalPopSize;
    nDied(ii,1) = sum(D.F(end,:)+D.Fv(end,:))*totalPopSize;
    maxHosp(ii,1) = max(sum(D.H' + D.Hv'))*totalPopSize;
    tMaxHosp(ii,1) = t(find(sum(D.H'+D.Hv')==max(sum(D.H'+D.Hv'))));
    
    % Save non-vaccinated results
    nInfs(ii,2) = sum(D.R(end,:))*totalPopSize;
    nHosps(ii,2) = sum(D.F(end,:)+D.Dis(end,:))*totalPopSize;
    nDied(ii,2) = sum(D.F(end,:))*totalPopSize;
    maxHosp(ii,2) = max(sum(D.H'))*totalPopSize;
    tMaxHosp(ii,2) = t(find(sum(D.H')==max(sum(D.H'))));
    
    % Save vaccinated results
    nInfs(ii,3) = sum(D.Rv(end,:))*totalPopSize;
    nHosps(ii,3) = sum(D.Fv(end,:)+D.Disv(end,:))*totalPopSize;
    nDied(ii,3) = sum(D.Fv(end,:))*totalPopSize;
    maxHosp(ii,3) = max(sum(D.Hv'))*totalPopSize;
    tMaxHosp(ii,3) = t(find(sum(D.Hv')==max(sum(D.Hv')), 1));
end

end