function Reff = getReff(P, V)

% getReff.m - calculate vaccinated reproduction number. See methods section 
% titled "Next generation matrix" in the main paper for details.
%
% Inputs:
%    P - Parameter structure obtained from getPar.m
%    V - 16xT matrix where V(i,t) is the proportion of age-group i that is
%    vaccinated at time t
%
% Outputs:
%    Reff - Tx1 vector where Reff(t) is the vaccinated reproduction number
%    at time t
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 26-12-2021


% Pre-allocate Reff vector
szV = size(V); Reff = zeros(szV(2), 1);

for ii = 1:szV(2)
    v = V(:,ii); % Extract current days vaccination vector (for simplicity)
    q = ((1 - P.VEi)*v) ./ (1 - P.VEi*v); % Calculate expected proportion of infections that are vaccinated
    NGM = P.tI * diag(1 - P.VEi * v) * diag(P.ui) * P.C' * diag(1 - P.VEt * q) * diag(P.IDR + P.tau*(1-P.IDR)); % Construct new un-scaled NGM
    Reff(ii) = max(abs(eig(P.u * NGM))); % Save current Reff
end


end