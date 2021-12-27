% resultsSEIR.m - reproduce SEIR model results
%
% Code to reproduce Figure 1, Table 2, Figure 2, Table 3, and Table 4
%
% Other m-files required: getReff.m, getVaccineRollout.m, getFinalEpidemicSize.m
% Subfunctions: getOpenBorderColumn()
%
% Author: Nic Steyn, Michael Plank
% Affiliation: Te Pūnaha Matatini
% Email: nicholas.steyn@auckland.ac.nz
% Last revision: 27-12-2021



%% Produce Figure 1 -------------------------------------------------------

clc
clear
close all
addpath("support")

% Inputs
r0_vals = [3.0, 4.5, 6.0];
nVax = (0:1e4:4.5e6)';
tMax = 730;
dailyImports = 5; 

% Get Results
Reff = zeros(length(nVax), 3, length(r0_vals));
nInfs = Reff; nHosps = Reff; nDied = Reff; maxHosp = Reff;

for rr = 1:length(r0_vals)

    % Setup
    r0 = r0_vals(rr);
    P = getPar(r0);
    [V, ageChanges] = getVaccineRollout(P, nVax, 0.9);
    m = dailyImports*P.popDist/sum(P.popDist);

    fprintf("Calculating results for R0 = %.1f...\n", r0)

    % Lower effectiveness
    P.VEi = 0.5; P.VEt = 0.4;
    Reff(:,1,rr) = getReff(P, V);
    [nInfsTemp, nHospsTemp, nDiedTemp, maxHospTemp, ~] = getFinalEpidemicSize(P, V, m, tMax);
    nInfs(:,1,rr) = nInfsTemp(:,1); nHosps(:,1,rr) = nHospsTemp(:,1);
    nDied(:,1,rr) = nDiedTemp(:,1); maxHosp(:,1,rr) = maxHospTemp(:,1);

    % Standard effectiveness
    P.VEi = 0.7; P.VEt = 0.5;
    Reff(:,2,rr) = getReff(P, V);
    [nInfsTemp, nHospsTemp, nDiedTemp, maxHospTemp, ~] = getFinalEpidemicSize(P, V, m, tMax);
    nInfs(:,2,rr) = nInfsTemp(:,1); nHosps(:,2,rr) = nHospsTemp(:,1);
    nDied(:,2,rr) = nDiedTemp(:,1); maxHosp(:,2,rr) = maxHospTemp(:,1);

    % Higher effectiveness
    P.VEi = 0.9; P.VEt = 0.5;
    Reff(:,3,rr) = getReff(P, V);
    [nInfsTemp, nHospsTemp, nDiedTemp, maxHospTemp, ~] = getFinalEpidemicSize(P, V, m, tMax);
    nInfs(:,3,rr) = nInfsTemp(:,1); nHosps(:,3,rr) = nHospsTemp(:,1);
    nDied(:,3,rr) = nDiedTemp(:,1); maxHosp(:,3,rr) = maxHospTemp(:,1);

end

fprintf("Calculating done..\n", r0)


% Plot results
col = [201 238 115; 74 189 133; 42 72 88]/255; % Colours for plotting
linewidth = 3;

titles = { ["Basic Reproduction Number: 3.0","(a)"], ["Basic Reproduction Number: 4.5", "(b)"], ["Basic Reproduction Number: 6.0", "(c)"],...
    "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)"};

for rr = 1:length(r0_vals)
    
    r0 = r0_vals(rr);

    % Plot Reff for r0 = r0_vals(rr)
    subplot(5, length(r0_vals), rr)

    p = plot(nVax/5e6, Reff(:,:,rr), 'LineWidth', linewidth);
    p(1).Color = col(1,:); p(2).Color = col(2,:); p(3).Color = col(3,:);
    ylim([0 r0]); xlim([0 0.9]);
    yline(1, 'HandleVisibility', 'off');
    title(titles{rr})

    xline(nVax(ageChanges(1)+1)/5e6,'HandleVisibility','off');
    h = text((nVax(ageChanges(1)+1)-5e4)/5e6, 0.2, "65+");
    set(h, 'Rotation', 90);
    xline(nVax(ageChanges(2)+1)/5e6,'HandleVisibility','off');
    h = text((nVax(ageChanges(2)+1)-5e4)/5e6, 0.2, "15+");
    set(h, 'Rotation', 90);

    if rr == 1; ylabel("Rv"); end

    % Plot infections
    subplot(5, length(r0_vals), rr + length(r0_vals))
    p = plot(nVax/5e6, nInfs(:,:,rr)/1e6, 'LineWidth', linewidth);
    p(1).Color = col(1,:); p(2).Color = col(2,:); p(3).Color = col(3,:);
    xlim([0 0.9])
    title(titles{rr + length(r0_vals)})

    xline(nVax(ageChanges(1)+1)/5e6,'HandleVisibility','off');
    xline(nVax(ageChanges(2)+1)/5e6,'HandleVisibility','off');

    if rr == 1;  ylabel(["Total Infections", "(millions)"]); end

    % Plot hospitalisations
    subplot(5, length(r0_vals), rr + 2*length(r0_vals))
    p = plot(nVax/5e6, nHosps(:,:,rr)/1e4, 'LineWidth', linewidth);
    p(1).Color = col(1,:); p(2).Color = col(2,:); p(3).Color = col(3,:);
    xlim([0 0.9])
    title(titles{rr + 2*length(r0_vals)})

    xline(nVax(ageChanges(1)+1)/5e6,'HandleVisibility','off');
    xline(nVax(ageChanges(2)+1)/5e6,'HandleVisibility','off');

    if rr == 1;  ylabel(["Total Hospitalisations", "(thousands)"]); end   

    % Plot deaths
    subplot(5, length(r0_vals), rr + 3*length(r0_vals))
    p = plot(nVax/5e6, nDied(:,:,rr)/1e4, 'LineWidth', linewidth);
    p(1).Color = col(1,:); p(2).Color = col(2,:); p(3).Color = col(3,:);
    xlim([0 0.9])
    title(titles{rr + 3*length(r0_vals)})

    xline(nVax(ageChanges(1)+1)/5e6,'HandleVisibility','off');
    xline(nVax(ageChanges(2)+1)/5e6,'HandleVisibility','off');

    if rr == 1;  ylabel(["Total Fatalities", "(thousands)"]); end  

    % Plot peak in hospital
    subplot(5, length(r0_vals), rr + 4*length(r0_vals))
    p = plot(nVax/5e6, maxHosp(:,:,rr)/1e4, 'LineWidth', linewidth);
    p(1).Color = col(1,:); p(2).Color = col(2,:); p(3).Color = col(3,:);
    xlim([0 0.9])
    title(titles{rr + 4*length(r0_vals)})

    xline(nVax(ageChanges(1)+1)/5e6,'HandleVisibility','off');
    xline(nVax(ageChanges(2)+1)/5e6,'HandleVisibility','off');

    if rr == 1;  ylabel(["Peak Hospital Occupancy", "(thousands)"]); end
    if rr == 2; xlabel("Proportion Fully Vaccinated"); end

end



%% Table 2 ----------------------------------------------------------------
% The table is produced entry by entry
% When required coverage is >90% user should find the value of maxVax such that maxVax =
% required vax coverage (at which point coverage will be equal in all age-groups)

clc
clear
close all
addpath("support")

% Inputs
r0 = 3.0;
vaxEff = 1; % = 1 for baseline, = 2 for lower, = 3 for higher
maxVax = 0.9;
VEi = [0.7 0.5 0.9]; VEt = [0.5 0.4 0.5];

% Setup
P = getPar(r0);
V = getVaccineRollout(P, (0:1e3:5e6)', maxVax);
P.VEi = VEi(vaxEff); P.VEt = VEt(vaxEff);

% Results
Reff = getReff(P, V);
thresh_ind = find(Reff<1, 1);
if ~isempty(thresh_ind)
    fprintf("%.2f%% of popn needs to be vaccinated\n", 100*sum(P.popDist.*V(:,thresh_ind)));
else
    fprintf("Cannot reach Reff < 1 with maxVax = %.2f and current vax eff.\n", maxVax)
end


%% Figure 2 ---------------------------------------------------------------

clc
clear
close all
addpath("support")

% Setup
r0_vec = 1:0.001:6;
V = getVaccineRollout(getPar(1), (0:1e3:4.5e6)', 0.9);
VEi = [0.7 0.5 0.9]; VEt = [0.5 0.4 0.5];
nVaxRequired = zeros(length(r0_vec), 3);

parfor ii = 1:length(r0_vec) % For each value of R0 (on x-axis)
    r0 = r0_vec(ii);
    P = getPar(r0);
    for jj = 1:3 % and for each vaccine effectiveness
        P.VEi = VEi(jj); P.VEt = VEt(jj);
        Reff = getReff(P, V);
        thresh_ind = find(Reff<1, 1); % find the point where Reff < 1
        if ~isempty(thresh_ind)
            nVaxRequired(ii, jj) = sum(P.popDist .* V(:,thresh_ind));
        else
            nVaxRequired(ii, jj) = inf;
        end
    end
end

% Plot
col = [201 238 115; 74 189 133; 42 72 88]/255;
P = getPar(1.0); % Just to get P.popDist

hold off
p = plot(r0_vec, nVaxRequired, 'LineWidth', 3);
p(1).Color = col(1,:); p(2).Color = col(2,:); p(3).Color = col(3,:);
hold on

yline(0.9*sum(P.popDist(4:end)), 'LineWidth', 1.5);
text(4.5, 0.75, "15-year-old threshold")

% Plot greatest R0 for all vax
xline(r0_vec(find(isinf(nVaxRequired(:,1)),1)), '--', 'Color', col(1,:), 'LineWidth', 1.5);
xline(r0_vec(find(isinf(nVaxRequired(:,2)),1)), '--', 'Color', col(2,:), 'LineWidth', 1.5);

% Plot greatest R0 for over-15s only
xline(r0_vec(find(nVaxRequired(:,1)>0.9*sum(P.popDist(4:end)),1)), ':', 'Color', col(1,:), 'LineWidth', 1.5);
xline(r0_vec(find(nVaxRequired(:,2)>0.9*sum(P.popDist(4:end)),1)), ':', 'Color', col(2,:), 'LineWidth', 1.5);
xline(r0_vec(find(nVaxRequired(:,3)>0.9*sum(P.popDist(4:end)),1)), ':', 'Color', col(3,:), 'LineWidth', 1.5);

ylim([0.4 0.9]);

ylabel("Proportion Required to be Vaccinated");
xlabel("Basic Reproduction Number");
legend(["Lower Effectiveness", "Baseline", "Higher Effectiveness","",""], 'Location', "SouthEast");


%% Table 3 & 4 ------------------------------------------------------------

clc
clear
close all
addpath('support')

% !!! SET VACCINE ALLOCATION HERE !!!
v = 0.9*ones(16,1); v(1:3) = 0; % For Table 3
% v = 0.9*ones(16,1); % For Table 4

% Setup
P = getPar(1.0); % Just to get P.popDist...
r0_vec = [3.0, 4.5, 6.0];
m = 5*P.popDist/sum(P.popDist); 
tMax = 730;

T = table;

for ii = 1:length(r0_vec)
    Ttemp = table;
    r0 = r0_vec(ii);
    P = getPar(r0);

    % Central effectiveness
    P.VEi = 0.7; P.VEt = 0.5;
    Ttemp.Central = getOpenBorderColumn(P, v, m, tMax);

    % Lower effectiveness
    P.VEi = 0.5; P.VEt = 0.4;
    Ttemp.Lower = getOpenBorderColumn(P, v, m, tMax);    
    
    % Higher effectiveness
    P.VEi = 0.9; P.VEt = 0.5;
    Ttemp.Higher = getOpenBorderColumn(P, v, m, tMax);

    Ttemp.Properties.RowNames = strcat(["Rv", "TotalInf", "Hosps", "Deaths", "PeakInHosp", "EndOf"], string(r0));
    
    T = [T; Ttemp];
end










%% Sub-functions ----------------------------------------------------------

function col = getOpenBorderColumn(P, v, m, tMax)
% Collates the columns for Tables 3 & 4

[nInfs, nHosps, nDied, maxHosp, tMaxHosp] = getFinalEpidemicSize(P, v, m, tMax);

col = strings(6, 1);
col(1) = sprintf("%.2f", getReff(P,v)); % Reff
col(2) = sprintf("%.2f (%.2f)", nInfs(1), nInfs(3)/nInfs(1));
col(3) = sprintf("%.2f (%.2f)", nHosps(1), nHosps(3)/nHosps(1));
col(4) = sprintf("%.2f (%.2f)", nDied(1), nDied(3)/nDied(1));
col(5) = sprintf("%.2f (after %d days)", [maxHosp(1), tMaxHosp(1)]);
col(6) = "";

end






