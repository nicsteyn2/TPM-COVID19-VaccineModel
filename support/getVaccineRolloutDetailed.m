function [V, ageChanges] = getVaccineRolloutDetailed(P, nVax, maxVax)

% getVaccineRolloutDetailed.m - like getVaccineRollout.m except vaccinates
% in finer age-groups (65+, 60-64, 55-59, ...., 15-19, 0-14).
%
% See getVaccineRollout.m for commented code
%
% Inputs:
%    P - Parameter structure obtained from getPar.m
%    nVax - Nx1 vector of total number of vaccine schedules
%    maxVax - value between 0 & 1 indicating maximum proportion of each
%    age-group that can be vaccinated
%
% Outputs:
%    V - 16xN matrix where V(i,j) is the proportion of age-group i
%    vaccinated when nVax(j) vaccine schedules are allocated
%    ageChanges - array of indices indicating when the next group is
%    vaccinated. This only works when nVax is monotonic. (Useful
%    for plotting)
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 26-12-2021

V = zeros(16, length(nVax));
group = zeros(length(nVax), 1);

groupSize = zeros(12,1);
groupSize(1) = sum(P.popCount(14:end));
groupSize(2:11) = P.popCount(13:-1:4);
groupSize(end) = sum(P.popCount(1:3));

for ii = 1:length(nVax)
    
    remainingVax = nVax(ii);
    nVaxByGroup = zeros(12,1);
    jj = 1;
    
    while (remainingVax > 0) && (jj <= 12)
        nVaxByGroup(jj) = min(remainingVax, maxVax*groupSize(jj));
        remainingVax = remainingVax - nVaxByGroup(jj);
        if remainingVax > 0
            jj = jj + 1;
        end
    end
    
    group(ii) = jj;
    
    V(14:end, ii) = nVaxByGroup(1)/groupSize(1);
    for kk = 4:13
        V(17-kk, ii) = nVaxByGroup(kk-2)/groupSize(kk-2);
    end
    V(1:3, ii) = nVaxByGroup(end)/groupSize(end);
    
end

ageChanges = find(diff(group)~=0);

end