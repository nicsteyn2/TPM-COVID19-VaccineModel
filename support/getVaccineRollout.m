function [V, ageChanges] = getVaccineRollout(P, nVax, maxVax)

% getVaccineRollout.m - calculate the age-structured vaccine rollout (65+,
% then 15-64, then 0-14) at the points described in vector nVax, with up to
% proportion maxVax in each age-group
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


V = zeros(16, length(nVax)); % Pre-allocate results matrix
group = zeros(length(nVax), 1); % Keep track of group that receives final vaccinations

% Calculate total size of each group
groupSize = zeros(3,1);
groupSize(1) = sum(P.popCount(14:end)); % Group 1: 65+
groupSize(2) = sum(P.popCount(4:13)); % Group 2: 15-64
groupSize(3) = sum(P.popCount(1:3)); % Group 3: 0-14

for ii = 1:length(nVax)
    
    remainingVax = nVax(ii); % Initially remainingVax is all vaccines at this step
    nVaxByGroup = zeros(3,1); % Save total number of vaccines allocated to each group

    jj = 1;
    while (remainingVax > 0) && (jj <= 3)
        nVaxByGroup(jj) = min(remainingVax, maxVax*groupSize(jj)); % Allocate either all vaccine to the group, or vaccines to everyone in group
        remainingVax = remainingVax - nVaxByGroup(jj); % Calculate remaining doses (only > 0 if there were more vaccines than people in group)
        if remainingVax > 0 % If there are leftovers, continue to next group
            jj = jj + 1;
        end
    end
    
    group(ii) = jj; % Save final group vaccinated
    
    % Convert from total counts to proportions
    V(14:end, ii) = nVaxByGroup(1)/groupSize(1);
    V(4:13, ii) = nVaxByGroup(2)/groupSize(2);
    V(1:3, ii) = nVaxByGroup(3)/groupSize(3);
    
end

ageChanges = find(diff(group)~=0); % Save when the final group vaccinated changed (assuming nVax is monotonic)

end