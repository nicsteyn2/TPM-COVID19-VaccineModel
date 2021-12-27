function [D, t] = runSEIR(P, ic, m, tMax)
% runSEIR.m - deterministic SEIR implementation
%
% Inputs:
%    P - Parameter structure obtained from getPar.m
%    m - 16x1 array where m(i) is the daily number of infected arrivals of age i
%    ic - 16x16 array of initial conditions (rows represent age-groups,
%    columns represent proportion of population in S, E, I, R compartments)
%    tMax - Maximum time of simulation in days
%
% Outputs:
%    D - Structure containing a tMax x 16 array for each compartment and 
%    three tMax x 16 arrays for cumulative cases, hosps, and fatalities.
%    t - a tMax x 1 array giving corresponding time
%
% Subfunctions: myODE
%
% Author: Nic Steyn, Michael Plank
% Te Pūnaha Matatini
% email: nicholas.steyn@auckland.ac.nz
% Last revision: 26-12-2021


%------------- Run Model --------------

opts = odeset('NonNegative', true, 'AbsTol', 1e-10); % Require all compartments to be non-negative
ic_reshaped = reshape(ic, [256, 1]); % Reshape initial conditions matrix into one row for ode45
[t, Y] = ode45(@(t, x)myODE(t, x, P, m), 0:tMax, ic_reshaped, opts); % Solve the ODE between time 0 and tMax at increments of 1 day


% ------------- Extract results --------------

D.S = Y(:,1:16); D.E = Y(:,17:32); D.I = Y(:,33:48); D.A = Y(:,49:64); D.R = Y(:,65:80);
D.Sv = Y(:,81:96); D.Ev = Y(:,97:112); D.Iv = Y(:,113:128); D.Av = Y(:,129:144); D.Rv = Y(:,145:160);
D.H = Y(:,161:176); D.Dis = Y(:, 177:192); D.F = Y(:,193:208);
D.Hv = Y(:,209:224); D.Disv = Y(:,225:240); D.Fv = Y(:,241:256);
D.All = Y;


end


%------------- ODE Function --------------

function f = myODE(t, x, P, m)

x = reshape(x, [16 16]); % Reshape state vector into human readable array
f = zeros(size(x)); % Output array

for ii = 1:16
    
    % Calculate the infection pressure acting on age group i
    LAM = (P.ui(ii)*P.u/P.popDist(ii)) *...
        sum( ((x(:,3) + P.tau * x(:,4)) + (1-P.VEt)*(x(:,8) + P.tau * x(:,9)) + ...
        P.tI * m./P.totalPopSize).*P.C(ii,:)' );
    
    % State transitions for non-vaccinated
    f(ii, 1) = -(LAM*x(ii,1)); % dS/dt
    f(ii, 2) = (LAM*x(ii,1)) - x(ii,2)/P.tE; % dE/dt
    f(ii, 3) = P.IDR(ii)*x(ii,2)/P.tE - x(ii,3)/P.tI; % dI/dt
    f(ii, 4) = (1-P.IDR(ii))*x(ii,2)/P.tE - x(ii,4)/P.tI; % dA/dt
    f(ii, 5) = (x(ii,3)+x(ii,4))/P.tI; % dR/dt
    
    % State transitions for vaccinated
    f(ii, 6) = -(LAM * x(ii,6)); % dSv/dt
    f(ii, 7) = (LAM * x(ii,6)) - x(ii,7)/P.tE; % dEv/dt
    f(ii, 8) = P.IDR(ii)*x(ii,7)/P.tE - x(ii,8)/P.tI; % dIv/dt
    f(ii, 9) = (1-P.IDR(ii))*x(ii,7)/P.tE - x(ii,9)/P.tI; % dAv/dt
    f(ii, 10) = (x(ii,8)+x(ii,9))/P.tI; % dRv/dt
    
    % State transitions for hospitalisation subroutine
    f(ii, 11) = (P.IHR(ii)*x(ii,2))/P.tE - x(ii,11)/P.tH; % dH/dt
    f(ii, 12) = (1 - P.IFR(ii)/P.IHR(ii))*x(ii,11)/P.tH; % dDischv/dt
    f(ii, 13) = (P.IFR(ii)/P.IHR(ii))*x(ii,11)/P.tH; % dF/dt
    f(ii, 14) = ((1-P.VEd)*P.IHR(ii)*x(ii,7))/P.tE - x(ii,14)/P.tH; % dHv/dt
    f(ii, 15) = (1 - P.IFR(ii)/P.IHR(ii))*x(ii,14)/P.tH; % dDischv/dt
    f(ii, 16) = (P.IFR(ii)/P.IHR(ii))*x(ii,14)/P.tH; % dF/dt

end

f = reshape(f, [256 1]); % Reshape into vector for ode45

end






