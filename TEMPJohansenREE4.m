clc; clear; close all
mrstVerbose on

%% --- MRST paths
restoredefaultpath; rehash toolboxcache
addpath(genpath('C:\Users\amira\Downloads\Advanced Subsurface Modelling\MATLAB\V8\JohansenCO2\mrst-2024b'))
rmpath(genpath(fullfile('C:\Users\amira\Downloads\Advanced Subsurface Modelling\MATLAB\V8\JohansenCO2\mrst-2024b\utils\octave_only')))
mrstModule add deckformat ad-core ad-props incomp mrst-gui

%% -------------------- Load Johansen Geometry --------------------
sector = 'Sector5';
fprintf(' Reading %s\n', sector);
grdecl = readGRDECL([sector '.grdecl']);

% Active cell mask (remove low-perm cells)
Kfile  = reshape(load([sector, '_Permeability.txt'])', [], 1);
pfile  = reshape(load([sector, '_Porosity.txt'])', [], 1);
grdecl.ACTNUM(Kfile < 0.11) = 0;

% Grid and geometry
G = processGRDECL(grdecl, 'Verbose', true);
G = computeGeometry(G);

%% -------------------- Rock & temperature fields --------------------
phi_raw       = pfile(G.cells.indexMap);
rock.poro     = max(0.01, min(0.35, phi_raw));
rock.perm     = 1.0 * bsxfun(@times, [1 1 0.5], Kfile(G.cells.indexMap)) * milli * darcy;
rock.rho_rock = 2650 * ones(G.cells.num, 1);

% Reservoir temperature (K) and reference
T_reservoir_C = 225 + 25 * exp(-0.1 * (G.cells.centroids(:,1)/1000).^2);
T_ref_C       = 25;
T_ref_K       = T_ref_C + 273.15;
T0_K          = T_reservoir_C + 273.15;

% Injection temperature
T_inj_C = 80;
T_inj_K = T_inj_C + 273.15;

% Clay content field (0..0.5)
clay_base  = 0.20; var_factor = 0.30;
rng(1234);
if isfield(G,'cartDims')
    try
        corr_len = max(G.cartDims)/10;
        raw = gaussianField(G.cartDims, [0, 1], [3, 3, 3], corr_len);
        v   = raw(G.cells.indexMap);
        v   = (v - mean(v(:))) / max(std(v(:)), 1e-12);
        rock.clay_content = clay_base * (1 + var_factor * v(:));
    catch
        rock.clay_content = clay_base * (1 + var_factor * randn(G.cells.num,1));
    end
else
    rock.clay_content = clay_base * (1 + var_factor * randn(G.cells.num,1));
end
rock.clay_content = min(max(rock.clay_content, 0), 0.5);

% Clay types and CEC (mol/kg rock), with pH effect
clay_props = struct( ...
    'illite',           struct('CEC', 0.100, 'fraction', 0.5), ...
    'kaolinite',        struct('CEC', 0.050, 'fraction', 0.3), ...
    'montmorillonite',  struct('CEC', 0.150, 'fraction', 0.2));

clay_type_idx = zeros(G.cells.num, 1);
if isfield(G,'cartDims')
    try
        tf = gaussianField(G.cartDims, [0, 1], [3,3,3], max(G.cartDims)/10);
        tf = tf(G.cells.indexMap);
        u  = (tf - min(tf)) ./ max(eps, (max(tf) - min(tf)));
        edges = [0, clay_props.illite.fraction, ...
                 clay_props.illite.fraction + clay_props.kaolinite.fraction, 1];
        clay_type_idx = discretize(u, edges);
    catch
        rv = rand(G.cells.num,1);
        clay_type_idx(rv <= clay_props.illite.fraction) = 1;
        clay_type_idx(rv > clay_props.illite.fraction & ...
                      rv <= clay_props.illite.fraction + clay_props.kaolinite.fraction) = 2;
        clay_type_idx(rv > clay_props.illite.fraction + clay_props.kaolinite.fraction) = 3;
    end
else
    rv = rand(G.cells.num,1);
    clay_type_idx(rv <= clay_props.illite.fraction) = 1;
    clay_type_idx(rv > clay_props.illite.fraction & ...
                  rv <= clay_props.illite.fraction + clay_props.kaolinite.fraction) = 2;
    clay_type_idx(rv > clay_props.illite.fraction + clay_props.kaolinite.fraction) = 3;
end
rock.clay_type = clay_type_idx;

rock.CEC = zeros(G.cells.num,1);
clay_names = fieldnames(clay_props);
for i = 1:G.cells.num
    tname = clay_names{clay_type_idx(i)};
    rock.CEC(i) = clay_props.(tname).CEC * rock.clay_content(i);
end

% pH attenuation of CEC
pH   = 5 * ones(G.cells.num,1);
pKa  = 6.5;
f_pH = 1 ./ (1 + 10.^(pH - pKa));
rock.CEC = rock.CEC .* f_pH;

% Base selectivities (for van 't Hoff later)
base_selectivity = struct('K_Li_Na', 2.0, 'K_Li_K', 5.0, 'K_Li_Mg', 8.0, 'K_Li_Ca', 10.0);

% Temperature dependence constants
R      = 8.314;    % J/(mol·K)
dH_ex  = -20e3;    % J/mol (equilibrium enthalpy for van 't Hoff)

% Transmissibilities
Trans = computeTrans(G, rock);

%% -------------------- Fluid --------------------
mu     = 0.3 * centi * poise;
rhoRef = 967 * kilogram/meter^3;
fluid  = initSingleFluid('mu', mu, 'rho', rhoRef);

%% -------------------- Initial Conditions --------------------
gravity on; g = gravity;
z  = G.cells.centroids(:,3);
z0 = max(z);
p0 = rhoRef * norm(g) * (z0 - z);
initState = initResSol(G, p0);

% Temperature as transported scalar (K)
initState.TK = T0_K;

% Aqueous concentrations (mol/kg water)
initState.c_Li = zeros(G.cells.num, 1);
initState.c_Na = 0.518  * ones(G.cells.num, 1);
initState.c_K  = 0.0036 * ones(G.cells.num, 1);
initState.c_Mg = 0.2416 * ones(G.cells.num, 1);
initState.c_Ca = 0.0076 * ones(G.cells.num, 1);

% Adsorbed equivalent fractions (sum to 1)
initState.N_Li = zeros(G.cells.num, 1);
initState.N_Na = 0.4 * ones(G.cells.num, 1);
initState.N_K  = 0.1 * ones(G.cells.num, 1);
initState.N_Mg = 0.3 * ones(G.cells.num, 1);
initState.N_Ca = 0.2 * ones(G.cells.num, 1);

%% -------------------- Wells (initial guess) --------------------
center_x = round(G.cartDims(1) / 2);
center_y = round(G.cartDims(2) / 2);

findWellCells = @(x_range, y_range, z_range) findValidWellCells(G, rock, x_range, y_range, z_range);

wc_inj_linear = findWellCells(center_x, center_y, 6:10);
if isempty(wc_inj_linear)
    validCells = find(rock.perm(:,1) > 1e-14 & ...
        z > min(z) + 0.3*(max(z)-min(z)) & z < min(z) + 0.7*(max(z)-min(z)));
    if ~isempty(validCells), wc_inj_linear = validCells(1:min(5, numel(validCells)));
    else, error('Cannot find suitable injection cells'); end
end

wc_prod_linear = findWellCells(center_x + 15, center_y, 6:10); % ensure some separation
if isempty(wc_prod_linear)
    validCells = find(rock.perm(:,1) > 1e-14 & ...
        z > min(z) + 0.3*(max(z)-min(z)) & z < min(z) + 0.7*(max(z)-min(z)));
    if ~isempty(validCells), wc_prod_linear = validCells(1:min(5, numel(validCells)));
    else, error('Cannot find suitable producer cells'); end
end

inj_rate  = 75/1000;  % m^3/s
prod_rate = 75/1000;  % m^3/s

W = [];
W = addWell(W, G, rock, wc_inj_linear,  'type', 'rate', 'val', inj_rate,  'comp_i', 1, 'Radius', 0.5*ft, 'name', 'Injector');
W = addWell(W, G, rock, wc_prod_linear, 'type', 'rate', 'val', -prod_rate, 'comp_i', 1, 'Radius', 0.5*ft, 'name', 'Producer');

fprintf('Initial injector cells: %s\n', mat2str(wc_inj_linear(:)'));
fprintf('Initial producer cells: %s\n', mat2str(wc_prod_linear(:)'));

bc = [];

%% ---------- Connectivity check & auto-relocate producer ----------
fprintf('\nPre-flight connectivity check (one pressure solve)...\n');
state0 = incompTPFA(initState, G, Trans, fluid, 'wells', W, 'bc', bc);

Nc    = G.cells.num;
Nn    = double(G.faces.neighbors);
intInx = all(Nn ~= 0, 2);
Nint  = Nn(intInx, :);
v0    = state0.flux(intInx);          % face flux (m^3/s), positive from Nint(:,1)->Nint(:,2)

% directed adjacency following flux direction
fpos = v0 > 0; fneg = v0 < 0;
rows = [Nint(fpos,1); Nint(fneg,2)];
cols = [Nint(fpos,2); Nint(fneg,1)];
A    = sparse(rows, cols, 1, Nc, Nc); % adjacency i->j if flux goes i->j

% BFS reachability from injector perforations
src   = unique(W(1).cells(:));
reach = false(Nc,1); queue = src(:); reach(queue) = true; head = 1;
while head <= numel(queue)
    i = queue(head); head = head + 1;
    nbrs = find(A(i,:));
    add  = nbrs(~reach(nbrs));
    reach(add) = true; queue = [queue; add(:)]; %#ok<AGROW>
end
tgt = unique(W(2).cells(:));
connected = any(reach(tgt));
fprintf('  Reachable producer from injector? %s\n', string(connected));

if ~connected
    fprintf('  Producer not reachable. Auto-relocating producer downstream...\n');
    cand = find(reach & rock.perm(:,1) > 1e-15);
    if isempty(cand), error('No reachable downstream permeable cells. Adjust ACTNUM or injector.'); end
    % rough graph distance
    depth = inf(Nc,1); depth(src) = 0; q = src(:); h=1;
    while h <= numel(q)
        i = q(h); h=h+1;
        nbrs = find(A(i,:));
        for j = nbrs(:)'
            if depth(j) > depth(i) + 1
                depth(j) = depth(i) + 1;
                q(end+1) = j; %#ok<AGROW>
            end
        end
    end
    [~,imax] = max(depth(cand)); farcell = cand(imax);
    c0 = G.cells.centroids(farcell,:); z0c = c0(3);
    dzOK = abs(G.cells.centroids(:,3) - z0c) < 0.25*(max(G.cells.centroids(:,3)) - min(G.cells.centroids(:,3)));
    d    = sqrt(sum((G.cells.centroids - c0).^2,2));
    pick = find(dzOK & reach & rock.perm(:,1) > 1e-15);
    [~,isrt] = sort(d(pick), 'ascend');
    wc_prod_new = pick(isrt(1:min(5,numel(isrt))));
    % rebuild producer
    W = W(1);
    W = addWell(W, G, rock, wc_prod_new, 'type','rate','val', -prod_rate, 'comp_i',1,'Radius',0.5*ft,'name','Producer');
    fprintf('  Producer moved. New cells: %s\n', mat2str(wc_prod_new(:)'));
    wc_prod_linear = wc_prod_new; % for later sampling
end
%% ---------- end connectivity check ----------

%% -------------------- Injection brine composition (mol/kg) --------------------
% Toggle scenarios here:
opt.injectLi = false;    % baseline: Li-free brine
opt.spike.Na = 1.0;      % multiply Na by this (e.g. 2, 5, ...)
opt.spike.K  = 1.0;
opt.spike.Mg = 1.0;
opt.spike.Ca = 1.0;

% background geothermal brine levels
Cin_bg = struct('Li', 0.0128, 'Na', 0.518, 'K', 0.0036, 'Mg', 0.2416, 'Ca', 0.146);

influxConc = struct();
influxConc.Li = Cin_bg.Li * double(opt.injectLi);
influxConc.Na = Cin_bg.Na * opt.spike.Na;
influxConc.K  = Cin_bg.K  * opt.spike.K;
influxConc.Mg = Cin_bg.Mg * opt.spike.Mg;
influxConc.Ca = Cin_bg.Ca * opt.spike.Ca;

fprintf('Injected brine (mol/kg): Li=%.4g Na=%.4g K=%.4g Mg=%.4g Ca=%.4g\n', ...
    influxConc.Li, influxConc.Na, influxConc.K, influxConc.Mg, influxConc.Ca);

%% -------------------- Boundary Conditions --------------------
% bc = [];

%% -------------------- Schedule --------------------
numSteps = 300;
totTime  = 3000 * day;
dt       = totTime / numSteps;
dt_min   = 0.5 * day;
dt_max   = 10  * day;
dt       = max(dt_min, min(dt, dt_max));
numSteps = ceil(totTime / dt);
fprintf('Using %d time steps of %.2f days each (nominal)\n', numSteps, convertTo(dt, day));

%% -------------------- State logs --------------------
states = cell(numSteps + 1,1);
states{1} = initState;

sol = struct('time', cell(numSteps + 1,1), 'pressure', [], ...
    'c_Li', [], 'c_Na', [], 'c_K', [], 'c_Mg', [], 'c_Ca', [], ...
    'N_Li', [], 'N_Na', [], 'N_K', [], 'N_Mg', [], 'N_Ca', []);
sol(1).time     = 0;
sol(1).pressure = initState.pressure / psia;
sol(1).c_Li     = initState.c_Li; sol(1).c_Na = initState.c_Na;
sol(1).c_K      = initState.c_K;  sol(1).c_Mg = initState.c_Mg;
sol(1).c_Ca     = initState.c_Ca;
sol(1).N_Li     = initState.N_Li; sol(1).N_Na = initState.N_Na;
sol(1).N_K      = initState.N_K;  sol(1).N_Mg = initState.N_Mg;
sol(1).N_Ca     = initState.N_Ca;

%% -------------------- Transport/Reaction parameters --------------------
k0_leach   = 1e-6;   % 1/s at T_ref (Arrhenius pre-exponential)
Ea_k       = 20e3;   % J/mol (activation energy)
cfl_factor = 0.6;    % for first-order upwind

Vp    = rock.poro .* G.cells.volumes;
rho_b = rock.rho_rock;

total_Li_dissolved  = zeros(numSteps + 1, 1);
total_Li_adsorbed   = zeros(numSteps + 1, 1);
recovery_efficiency = zeros(numSteps + 1, 1);

rhoRef = convertFrom(rhoRef, kilogram/meter^3); % keep numeric for ops

total_Li_dissolved(1) = sum(rhoRef * states{1}.c_Li .* Vp);
total_Li_adsorbed(1)  = sum(rho_b .* (1 - rock.poro) .* rock.CEC .* states{1}.N_Li .* G.cells.volumes);

t_accum = 0;

%% -------------------- Discrete operators & diffusion geometry --------------------
Nn     = double(G.faces.neighbors);
intInx = all(Nn ~= 0, 2);
Nint   = Nn(intInx, :);
nf     = size(Nint,1);
Cmat   = sparse([(1:nf)';(1:nf)'], Nint, ones(nf,1)*[-1 1], nf, G.cells.num);

Af  = G.faces.areas(intInx);
i1  = Nint(:,1); i2 = Nint(:,2);
d12 = sqrt(sum((G.cells.centroids(i2,:) - G.cells.centroids(i1,:)).^2, 2));
phi_face = 0.5*(rock.poro(i1) + rock.poro(i2));

% Solute diffusion/dispersion parameters
D_m    = 1e-9;     % molecular diffusion [m^2/s]
alphaL = 20.0;     % longitudinal dispersivity [m]

% Thermal diffusion/dispersion parameters
D_T_m    = 1.0e-6; % thermal diffusivity [m^2/s]
alphaL_T = 5.0;    % thermal dispersivity [m]

%% -------------------- Main Time Loop --------------------
cations = {'Li','Na','K','Mg','Ca'};
zval    = [1, 1, 1, 2, 2];
Nc      = G.cells.num;

for step = 1:numSteps
    fprintf('\nTime step %d\n', step);

    % Pressure solve
    state = incompTPFA(states{step}, G, Trans, fluid, 'wells', W, 'bc', bc); %#ok<NBRAK>
    % bracket fix for MATLAB: use parentheses, not brackets:
    % state = incompTPFA(states{step}, G, Trans, fluid, 'wells', W, 'bc', bc);

    % Carry scalar fields
    for j = 1:numel(cations)
        nm = cations{j};
        state.(['c_',nm]) = states{step}.(['c_',nm]);
        state.(['N_',nm]) = states{step}.(['N_',nm]);
    end
    state.TK = states{step}.TK;

    % Well perforation fluxes
    q_inj  = zeros(Nc, 1); q_inj(W(1).cells)  = state.wellSol(1).flux;
    q_prod = zeros(Nc, 1); q_prod(W(2).cells) = state.wellSol(2).flux;
    fprintf('  Well balance (should ~0): %.3e m^3/s\n', sum(q_inj) + sum(q_prod));

    % Face fluxes and solute diffusion weights
    v = state.flux(intInx);
    u_face = abs(v) ./ Af;
    D_eff  = D_m + alphaL * u_face;
    Diff_w = phi_face .* D_eff .* Af ./ max(d12, 1e-12);

    % Reaction limiter using current temperature
    TempK_prev   = state.TK;
    k_leach_prev = k0_leach .* exp(-Ea_k./R .* (1./TempK_prev - 1./T_ref_K));
    dt_reac      = 0.2 ./ max(k_leach_prev, 1e-18);

    % Advection CFL and diffusion limiters
    F_out   = full(max(0, -Cmat)'*min(v,0) + max(0, Cmat)'*max(v,0));
    F_out   = F_out + abs(q_prod);
    dt_cfl  = cfl_factor * Vp ./ max(F_out, 1e-18);
    dt_diff = 0.45 * min(Vp) / max(Diff_w, [], 'omitnan');

    dt_actual = min([dt, max(dt_min, min(dt_max, min([dt_cfl; dt_diff; dt_reac])) )]);

    % ---------- Minimal adaptive retry loop ----------
    dt_try     = dt_actual;
    base_state = state;
    accepted   = false;
    mbe = NaN;

    for retry = 1:5
        state = base_state;

        % Temperature advection + diffusion
        TempK = state.TK;
        flag  = v > 0;

        upT    = TempK(Nint(:,1)).*double(flag) + TempK(Nint(:,2)).*double(~flag);
        F_T    = upT .* v;
        divF_T = -Cmat' * F_T;

        u_face_T  = abs(v) ./ Af;
        D_eff_T   = D_T_m + alphaL_T * u_face_T;
        Diff_w_T  = phi_face .* D_eff_T .* Af ./ max(d12, 1e-12);
        gradT     = Cmat * TempK;
        Fdiff_T   = - Diff_w_T .* gradT;
        divDiffT  = - Cmat' * Fdiff_T;

        S_T_inj  = q_inj * T_inj_K;
        S_T_prod = q_prod .* TempK;
        rhsT     = -divF_T + S_T_inj + S_T_prod + divDiffT;
        TempK    = TempK + dt_try .* (rhsT ./ max(Vp, 1e-30));
        TempK    = max(250, min(500, TempK));

        % T-dependent selectivities (van 't Hoff) and kinetics (Arrhenius)
        temp_fac_cell = exp(-dH_ex./R .* (1./TempK - 1./T_ref_K));
        K_Li_Na_T = base_selectivity.K_Li_Na .* temp_fac_cell;
        K_Li_K_T  = base_selectivity.K_Li_K  .* temp_fac_cell;
        K_Li_Mg_T = base_selectivity.K_Li_Mg .* temp_fac_cell;
        K_Li_Ca_T = base_selectivity.K_Li_Ca .* temp_fac_cell;
        Kvec_T = [K_Li_Na_T(:), K_Li_K_T(:), K_Li_Mg_T(:), K_Li_Ca_T(:)];
        k_leach_cell = k0_leach .* exp(-Ea_k./R .* (1./TempK - 1./T_ref_K));

        % Equilibrium and exact kinetic relaxation
        cMat  = [state.c_Li, state.c_Na, state.c_K, state.c_Mg, state.c_Ca];
        f_prev = [state.N_Li, state.N_Na, state.N_K, state.N_Mg, state.N_Ca];

        f_eq = zeros(Nc, 5);
        for iCell = 1:Nc
            f_eq(iCell,:) = solveGTEqFractions_Fixed(cMat(iCell,:), Kvec_T(iCell,:), zval, f_prev(iCell,:));
        end
        theta   = exp(-k_leach_cell .* dt_try);
        fMatNew = zeros(Nc, 5);
        for jj = 1:5
            fMatNew(:,jj) = f_prev(:,jj).*theta + f_eq(:,jj).*(1 - theta);
        end
        fMatNew = renormaliseFractions(fMatNew);

        % Mass moved to solid (mol per cell per species)
        Delta_ads_mol = rock.rho_rock .* (1 - rock.poro) .* rock.CEC .* (fMatNew - f_prev) .* G.cells.volumes;

        % Solute transport for each cation
        neg_flag = false;
        for j = 1:5
            nm = cations{j};
            c  = state.(['c_',nm]);

            upc    = c(Nint(:,1)).*double(flag) + c(Nint(:,2)).*double(~flag);
            F_adv  = rhoRef * upc .* v;
            divAdv = -Cmat' * F_adv;

            gradc   = Cmat * c;
            Fdiff   = - rhoRef .* Diff_w .* gradc;
            divDiff = - Cmat' * Fdiff;

            S_inj  = rhoRef .* q_inj * influxConc.(nm);
            S_prod = rhoRef .* q_prod .* c;
            S_well = S_inj + S_prod;

            S_ads  = - Delta_ads_mol(:, j) ./ max(dt_try, eps);

            rhs   = -divAdv + S_well + S_ads + divDiff;
            denom = rhoRef * Vp;
            c_new = c + dt_try .* (rhs ./ max(denom, 1e-30));
            if any(c_new < -1e-10), neg_flag = true; end
            state.(['c_',nm]) = max(c_new, 0);
        end

        % Commit solid and temperature
        state.N_Li = fMatNew(:,1); state.N_Na = fMatNew(:,2); state.N_K  = fMatNew(:,3);
        state.N_Mg = fMatNew(:,4); state.N_Ca = fMatNew(:,5); state.TK   = TempK;

        % Per-substep Li mass balance
        tin   = sum(q_inj)  * rhoRef * influxConc.Li * dt_try;
        tout  = abs(sum(q_prod)) * rhoRef * mean(state.c_Li(wc_prod_linear)) * dt_try;
        diss0 = sum(rhoRef * base_state.c_Li .* Vp);
        ads0  = sum(rock.rho_rock .* (1 - rock.poro) .* rock.CEC .* base_state.N_Li .* G.cells.volumes);
        diss1 = sum(rhoRef * state.c_Li .* Vp);
        ads1  = sum(rock.rho_rock .* (1 - rock.poro) .* rock.CEC .* state.N_Li .* G.cells.volumes);
        mbe   = abs(tin - tout - (diss1 + ads1 - diss0 - ads0)) / max(tin + 1e-30, 1);

        if ~neg_flag && (mbe <= 0.02 || tin == 0)
            accepted = true; break
        else
            dt_try = 0.5 * dt_try;
        end
    end

    if ~accepted
        warning('Accepted step with relaxed checks after retries. mbe=%.3f', mbe);
    end

    % Advance time by the accepted dt
    t_accum = t_accum + dt_try;

    %% Diagnostics and logs
    total_mass_in  = sum(q_inj)  * rhoRef * influxConc.Li * dt_try;
    prod_c_Li      = mean(state.c_Li(wc_prod_linear));
    total_mass_out = abs(sum(q_prod)) * rhoRef * prod_c_Li * dt_try;

    dissolved_Li      = sum(rhoRef * state.c_Li .* Vp);
    adsorbed_Li       = sum(rock.rho_rock .* (1 - rock.poro) .* rock.CEC .* state.N_Li .* G.cells.volumes);
    dissolved_Li_prev = sum(rhoRef * states{step}.c_Li .* Vp);
    adsorbed_Li_prev  = sum(rock.rho_rock .* (1 - rock.poro) .* rock.CEC .* states{step}.N_Li .* G.cells.volumes);
    mass_stored       = (dissolved_Li + adsorbed_Li) - (dissolved_Li_prev + adsorbed_Li_prev);

    fprintf('  Li mass: In=%.2e, Out=%.2e, Stored=%.2e (mol)\n', total_mass_in, total_mass_out, mass_stored);
    fprintf('  Producer c_Li: %.4e mol/kg\n', prod_c_Li);

    total_Li_dissolved(step + 1)  = dissolved_Li;
    total_Li_adsorbed(step + 1)   = adsorbed_Li;
    recovery_efficiency(step + 1) = total_mass_out / max(total_mass_in, eps) * 100;

    states{step + 1}       = state;
    sol(step + 1).time     = t_accum;
    sol(step + 1).pressure = state.pressure / psia;
    sol(step + 1).c_Li     = state.c_Li;  sol(step + 1).c_Na = state.c_Na;
    sol(step + 1).c_K      = state.c_K;   sol(step + 1).c_Mg = state.c_Mg;
    sol(step + 1).c_Ca     = state.c_Ca;
    sol(step + 1).N_Li     = state.N_Li;  sol(step + 1).N_Na = state.N_Na;
    sol(step + 1).N_K      = state.N_K;   sol(step + 1).N_Mg = state.N_Mg;
    sol(step + 1).N_Ca     = state.N_Ca;
end

%% -------------------- Analysis --------------------
prodConc     = zeros(numSteps + 1, 1);
prodAdsorbed = zeros(numSteps + 1, 1);
for i = 1:numSteps + 1
    prodConc(i)     = mean(states{i}.c_Li(wc_prod_linear));
    prodAdsorbed(i) = mean(states{i}.N_Li(wc_prod_linear) .* rock.CEC(wc_prod_linear));
end

breakthroughThreshold = 0.01 * max(influxConc.Li, eps);
idx = find(prodConc >= breakthroughThreshold, 1);
if isempty(idx)
    breakthroughTime = NaN; fprintf('Breakthrough not reached within simulated time.\n');
else
    breakthroughTime = sol(idx).time;
    fprintf('Breakthrough time: %.2f days\n', convertTo(breakthroughTime, day));
end
final_recovery = recovery_efficiency(end);

%% -------------------- Visualisation --------------------
t_days = arrayfun(@(s) convertTo(s.time, day), sol);

figure('Position', [100, 100, 1200, 800]);
subplot(2,3,1);
plot(t_days, prodConc, 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Li Conc (mol/kg)'); title('Dissolved Li at Producer'); grid on;
if ~isnan(breakthroughTime), hold on; plot(convertTo(breakthroughTime, day), breakthroughThreshold, 'ro'); end

subplot(2,3,2);
plot(t_days, prodAdsorbed, 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Adsorbed Li (mol/kg rock)'); title('Adsorbed Li near Producer'); grid on;

subplot(2,3,3);
plot(t_days, total_Li_dissolved, 'LineWidth', 2); hold on;
plot(t_days, total_Li_adsorbed, 'LineWidth', 2);
plot(t_days, total_Li_dissolved + total_Li_adsorbed, 'k--');
xlabel('Time (days)'); ylabel('Li Mass (mol)'); title('System Li Mass Balance'); legend('Dissolved','Adsorbed','Total'); grid on;

subplot(2,3,4);
cLi = max(states{end}.c_Li, 1e-12);
plotCellData(G, log10(cLi), 'EdgeAlpha', 0.1); plotWell(G, W);
title('Final log_{10}(Li)'); colorbar;
caxis([log10(1e-6), log10(max(max(cLi), max(influxConc.Li,1e-6)))]); view(-63, 68); axis tight off;

subplot(2,3,5);
plotCellData(G, states{end}.N_Li .* rock.CEC, 'EdgeAlpha', 0.1); plotWell(G, W);
title('Final Adsorbed Li (mol/kg rock)'); colorbar; view(-63, 68); axis tight off;

subplot(2,3,6);
plotCellData(G, states{end}.TK - 273.15, 'EdgeAlpha', 0.1); plotWell(G, W);
title('Temperature (°C)'); colorbar; view(-63, 68); axis tight off;

try
    mrstModule add mrst-gui
    figure; plotToolbar(G, states, 'field', 'c_Li'); view(-63, 68); colorbar;
    title('Dissolved Li Evolution');
catch
end

%% -------------------- Summary --------------------
fprintf('\n=== REACTIVE TRANSPORT SIMULATION SUMMARY ===\n');
fprintf('Total simulated time: %.2f days\n', convertTo(sol(end).time, day));
fprintf('Number of steps: %d\n', numSteps);
fprintf('Grid: %d x %d x %d | Active cells: %d\n', G.cartDims(1), G.cartDims(2), G.cartDims(3), G.cells.num);
fprintf('Injection rate: %.1f stb/day | Production rate: %.1f stb/day\n', ...
    convertTo(abs(W(1).val), stb/day), convertTo(abs(W(2).val), stb/day));
fprintf('Injected brine (mol/kg): Li=%.4g Na=%.4g K=%.4g Mg=%.4g Ca=%.4g\n', ...
    influxConc.Li, influxConc.Na, influxConc.K, influxConc.Mg, influxConc.Ca);
fprintf('Final producer Li concentration: %.4e mol/kg\n', prodConc(end));
fprintf('Average clay content: %.2f %% | Average CEC: %.3f mol/kg rock\n', mean(rock.clay_content)*100, mean(rock.CEC));
if ~isnan(breakthroughTime), fprintf('Breakthrough time: %.2f days\n', convertTo(breakthroughTime, day));
else, fprintf('Breakthrough time: Not reached\n'); end
fprintf('Final recovery efficiency (Li): %.2f %%\n', final_recovery);
fprintf('============================================\n');

%% -------------------- Helper Functions --------------------
function cells = findValidWellCells(G, rock, x_range, y_range, z_range)
    cells = [];
    if isfield(G,'cartDims')
        try
            wc = false(G.cartDims);
            wc(x_range, y_range, z_range) = true;
            candidate_cells = find(wc(G.cells.indexMap));
            perm_threshold  = 1e-14; % ~10 mD
            valid_idx = rock.perm(candidate_cells, 1) > perm_threshold;
            cells = candidate_cells(valid_idx);
            if isempty(cells), cells = candidate_cells; end
        catch
        end
    end
end

function f_eq = solveGTEqFractions_Fixed(c, K_Na, z, f_guess)
    % c = [c_Li c_Na c_K c_Mg c_Ca], z = [1 1 1 2 2], Na is reference
    Na_idx = 2; c = max(c, 1e-12);
    K_Li_Na = K_Na(1); K_Li_K = K_Na(2); K_Li_Mg = K_Na(3); K_Li_Ca = K_Na(4);
    K_vs_Na = [K_Li_Na, 1.0, K_Li_Na / K_Li_K, K_Li_Na / K_Li_Mg, K_Li_Na / K_Li_Ca];

    alpha = zeros(1, 5);
    for i = 1:5
        if i == Na_idx, continue; end
        alpha(i) = K_vs_Na(i) * (c(i)^z(Na_idx)) / (c(Na_idx)^z(i));
    end
    % Quadratic for f_Na: a*f_Na^2 + b*f_Na - 1 = 0
    a = alpha(4)+alpha(5);
    b = alpha(1)+1+alpha(3);
    rootsNa = roots([a, b, -1]);
    f_Na = rootsNa(rootsNa > 0 & rootsNa < 1);
    if isempty(f_Na)
        f_eq = f_guess / max(sum(f_guess), 1e-12);
        return
    end
    f_Na = f_Na(1);

    f_eq = zeros(1,5);
    f_eq(Na_idx) = f_Na;
    f_eq(1) = alpha(1) * f_Na^(z(1)/z(Na_idx)); % Li
    f_eq(3) = alpha(3) * f_Na^(z(3)/z(Na_idx)); % K
    f_eq(4) = alpha(4) * f_Na^(z(4)/z(Na_idx)); % Mg
    f_eq(5) = alpha(5) * f_Na^(z(5)/z(Na_idx)); % Ca

    f_eq = max(0, f_eq); f_eq = f_eq / sum(f_eq);
end

function F = renormaliseFractions(F)
    s = sum(F, 2); s = max(s, 1e-12);
    F = F ./ s; F = max(0, min(1, F));
end
