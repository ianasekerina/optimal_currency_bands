%% =========================================================================
%  OptimalBands.m
%  Replication code for Thesis: "Optimal Currency Bands" 
%
%  This script reproduces all figures and calibration results in the thesis.
%  The model builds on Krugman (1991): a monetary authority chooses a
%  symmetric exchange-rate band [−b, b] to minimise steady-state exchange-
%  rate variance plus proportional FX-intervention costs.
%
%  HOW TO USE
%  ----------
%  1. Set DATA_DIR (line ~48) to the folder containing the two ECB CSV files.
%  2. Run the script.  All figures are produced; calibration results are
%     printed to the console and saved as structured output variables.
%
%  OUTPUT VARIABLES (available in the workspace after running)
%  -----------------------------------------------------------
%  res_2pct   — calibration results for the ±2.25 % ERM band
%  res_6pct   — calibration results for the ±6 %   ERM band
%  cf1, cf2   — counterfactual forward-solve results
%
%  DEPENDENCIES
%  ------------
%  MATLAB R2019b or later (uses `tiledlayout`, `xline`, `yline`).
%  No additional toolboxes required; all root-finding uses built-in `fzero`.
%
%  REFERENCES
%  ----------
%  Krugman, P. (1991). "Target Zones and Exchange Rate Dynamics."
%    QJE, 106(3), 669–682.
%  Miller, M. & Zhang, L. (1996). "Optimal target zones."
%    JEDC, 20(9), 1641–1660.
%  Coenen, G. & Vega, J.-L. (2001). "The Demand for M3 in the Euro Area."
%    JAE, 16(6), 727–748.
%
%  Author: Iana Sekerina, 2026.
%  AI tools were used for code debugging and editorial assistance only.
% =========================================================================

clear; close all;

% -------------------------------------------------------------------------
%  DATA DIRECTORY — UPDATE THIS PATH BEFORE RUNNING
%  The folder should contain:
%    ert_h_eur_d__custom_19839545_linear.csv  (ECB daily ITL/ECU, main)
%    ert_h_eur_d__custom_19839614_linear.csv  (ECB daily ITL/ECU, supplement)
% -------------------------------------------------------------------------
DATA_DIR = '/Users/yanasekerina/Downloads/';

% =========================================================================
%  SECTION 0 — COLOUR PALETTE
%  A consistent, colour-blind-friendly palette used throughout all figures.
% =========================================================================

COL.obs    = [0,   0,   0  ] / 255;        % black   — observed exchange rate
COL.band6  = [230, 159,  0 ] / 255;        % orange  — ±6 % ERM band (wide
COL.band2  = [  0, 114, 178] / 255;        % blue    — ±2.25 % ERM band (narrow)
COL.band15 = [  0, 158, 115] / 255;        % teal    — ±15 % / generic "large"
COL.parity = [0.35, 0.35, 0.35];           % grey    — central parities / baselines

% Aliases used in the approximation-error figures (Section 2)
clr_approx_small = COL.band15;   % narrow-band approximation curves
clr_approx_large = COL.band6;    % wide-band  approximation curves
clr_baseline     = COL.parity;   % baseline / zero-error reference lines

% =========================================================================
%  SECTION 1 — KRUGMAN EXCHANGE RATE FUNCTION  (Figure 1 in thesis)
%
%  Proposition 1: under a credible symmetric band [−b, b] the exchange rate
%  is the S-shaped function
%
%      f(x) = (1/β) [ x − (1/η) sinh(ηx)/cosh(ηX) ],   η = sqrt(2β/σ²)
%
%  where X (fundamental boundary) is linked to b by:
%      b = (1/ηβ)(ηX − tanh(ηX))
%
%  The S-curve lies strictly inside the free-float line s = x/β.
% =========================================================================

fprintf('\n=== SECTION 1 : Krugman exchange rate function ===\n');

% Model parameters for the illustration 
beta_fig1 = 1.0;     % inverse of money-demand semi-elasticity (α = 1)
sigma_fig1 = 0.1;    % fundamental volatility σ
X_fig1 = 0.15;       % fundamental boundary (half-width of x's domain)

eta_fig1 = sqrt(2 * beta_fig1 / sigma_fig1^2);   % curvature parameter η

% Exchange rate function (eq. 5) and free-float benchmark 
f_tz    = @(x)  (1/beta_fig1) .* (x - (1/eta_fig1) .* ...
                sinh(eta_fig1 .* x) ./ cosh(eta_fig1 .* X_fig1));
f_float = @(x)  x ./ beta_fig1;

% Band half-width implied by X (eq. 6) 
b_fig1 = (1 / (eta_fig1 * beta_fig1)) * (eta_fig1 * X_fig1 - tanh(eta_fig1 * X_fig1));
fprintf('Illustration parameters: β=%.2f, σ=%.2f, X=%.2f  –  b=%.4f\n', ...
    beta_fig1, sigma_fig1, X_fig1, b_fig1);

x_grid = linspace(-X_fig1, X_fig1, 1000);

figure('Position', [100 100 1100 900], 'Color', 'w');
plot(x_grid, f_tz(x_grid),    'LineWidth', 3, 'Color', COL.band15); hold on;
plot(x_grid, f_float(x_grid), 'LineWidth', 3, 'Color', COL.band2);
yline( b_fig1, 'k', 'LineWidth', 2);
yline(-b_fig1, 'k', 'LineWidth', 2);
xline( X_fig1, 'LineWidth', 2, 'Color', COL.band15);
xline(-X_fig1, 'LineWidth', 2, 'Color', COL.band15);
xlim([-X_fig1 - 0.01,  X_fig1 + 0.01]);
ylim([-b_fig1 - 0.001,  b_fig1 + 0.001]);
xlabel('Fundamentals $x$',   'Interpreter', 'latex', 'FontSize', 20);
ylabel('Exchange rate $s$',  'Interpreter', 'latex', 'FontSize', 20);
legend({'$f(x)$', '$x/\beta$', '$\pm b$', '$\pm X$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 20);
set(gca, 'FontSize', 16);
grid on;

% =========================================================================
%  SECTION 2 — APPROXIMATION ACCURACY  (Figure 2 & Appendix Figure 8)
%
%  We compare the exact optimal band X* (from FOC eq. 16) with its
%  narrow-band (eq. 17) and wide-band (eq. 18) approximations.
%
%  Two sub-figures are produced:
%   (a) Log deviation as a function of θ = κ/σ         – Figure 2
%   (b) Log deviation as a function of κ and σ separately – Appendix Fig 8
% =========================================================================

fprintf('\n=== SECTION 2 : Approximation accuracy ===\n');

% Baseline β for this section (Coenen & Vega 2001 money-demand estimate)
alpha_base = 0.87;
beta  = 1 / alpha_base;   % β = 1/α ≈ 1.1494;  used for all remaining sections

% Homogeneous g-function (Corollary 2): g(v) = κ/σ at the optimum 
%  v := X/σ;  g strictly increasing, so the equation g(v) = θ has a unique
%  solution for each θ = κ/σ > 0.
g_fun = @(v) ...
    (4 * v.^3) ./ (3 * beta^2) ...
    - 5 * tanh(sqrt(2*beta) .* v) ./ (2 * sqrt(2) * beta^(7/2)) ...
    + 5 * v .* sech(sqrt(2*beta) .* v).^2 ./ (2 * beta^3) ...
    + sqrt(2) .* v.^2 .* tanh(sqrt(2*beta) .* v) ...
        .* sech(sqrt(2*beta) .* v).^2 ./ beta^(5/2);

% Approximate v* as a function of θ (from Proposition 3) 
phi_narrow = @(theta) (105/272)^(1/7) .* theta.^(1/7);   % narrow-band
phi_wide   = @(theta) (3 * beta^2 / 4)^(1/3) .* theta.^(1/3); % wide-band

% (a) Figure 2: log deviation vs θ = κ/σ 
%  High θ  –  narrow-band regime (high cost or low volatility)
%  Low  θ  –  wide-band regime   (low cost or high volatility)
theta_narrow_grid = linspace(0.000001,  0.01,  10000);
theta_wide_grid   = linspace(4, 100,  10000);

v_exact_narrow = nan(size(theta_narrow_grid));
v_exact_wide   = nan(size(theta_wide_grid));

for i = 1:length(theta_narrow_grid)
    th = theta_narrow_grid(i);
    try; v_exact_narrow(i) = fzero(@(v) g_fun(v) - th, [1e-7, 100]); catch; end
end
for i = 1:length(theta_wide_grid)
    th = theta_wide_grid(i);
    try; v_exact_wide(i) = fzero(@(v) g_fun(v) - th, [1e-7, 1000]); catch; end
end

logdev_narrow_theta = 100 * log(phi_narrow(theta_narrow_grid) ./ v_exact_narrow);
logdev_wide_theta   = 100 * log(phi_wide(theta_wide_grid)     ./ v_exact_wide);


fig2 = figure('Position', [100 100 1300 520], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

nexttile;
hold on;
plot(theta_narrow_grid, logdev_narrow_theta, 'LineWidth', 2.5, 'Color', clr_approx_small);
yline(0, '--', 'Color', clr_baseline, 'LineWidth', 1.5);
hold off;
xlabel('$\theta$',  'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex', 'FontSize', 14);
title('Narrow band approximation', 'Interpreter', 'latex', 'FontSize', 20);
xlim([min(theta_narrow_grid), max(theta_narrow_grid)]);
ylim([-10, 10]);
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 16);
grid on; grid minor;

nexttile;
hold on;
plot(theta_wide_grid, logdev_wide_theta, 'LineWidth', 2.5, 'Color', clr_approx_large);
yline(0, '--', 'Color', clr_baseline, 'LineWidth', 1.5);
xline(theta_calibrated, ':', 'Color', 'black', 'LineWidth', 2, ...
    'Label', 'ERM I calibration ($\theta\approx0.47$)', 'Interpreter', 'latex', ...
    'LabelOrientation', 'horizontal', 'FontSize', 13, 'HandleVisibility', 'off');
hold off;
xlabel('$\theta$',  'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex', 'FontSize', 14);
title('Wide band approximation', 'Interpreter', 'latex', 'FontSize', 20);
xlim([min(theta_wide_grid), max(theta_wide_grid)]);
ylim([-10, 10]);
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 16);
grid on; grid minor;

% Combined Figure 2: both approximations, log-theta axis, calibration marked 
theta_grid = logspace(-4, log10(200), 10000);
v_exact_grid = nan(size(theta_grid));
for i = 1:length(theta_grid)
    try; v_exact_grid(i) = fzero(@(v) g_fun(v) - theta_grid(i), [1e-8, 1000]); catch; end
end
logdev_narrow = 100 * log(phi_narrow(theta_grid) ./ v_exact_grid);
logdev_wide   = 100 * log(phi_wide(theta_grid)   ./ v_exact_grid);

theta_calibrated = 0.47;  % both ERM I regimes land here (0.4720, 0.4717)

fig2 = figure('Position', [100 100 900 550], 'Color', 'w');
hold on;
plot(theta_grid, logdev_narrow, 'LineWidth', 2.5, 'Color', clr_approx_small, ...
    'DisplayName', 'Narrow-band approximation');
plot(theta_grid, logdev_wide,   'LineWidth', 2.5, 'Color', clr_approx_large, ...
    'DisplayName', 'Wide-band approximation');
yline(0, '--', 'Color', clr_baseline, 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(theta_calibrated, ':', 'Color', 'black', 'LineWidth', 2, ...
    'Label', 'ERM I calibration ($\theta\approx0.47$)', 'Interpreter', 'latex', ...
    'LabelOrientation', 'horizontal', 'FontSize', 13, 'HandleVisibility', 'off');
hold off;
set(gca, 'XScale', 'log', 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 13);
ylim([-30, 10]);
grid on; grid minor;
% (b) Appendix Figure 8: log deviation varying κ and σ separately 
%  Parameters held fixed while the other varies:
kappa_fix_narrow = 0.0001;   sigma_fix_narrow = 0.045; % narrow-band reference
kappa_fix_wide   = 5.000;   sigma_fix_wide   = 0.300; % wide-band  reference

kappa_grid_narrow = linspace(0.6*kappa_fix_narrow, 1.4*kappa_fix_narrow, 1000);  
sigma_grid_narrow = linspace(0.030, 0.060, 1000);                                
kappa_grid_wide   = linspace(0.6*kappa_fix_wide,   1.4*kappa_fix_wide,   1000);   
sigma_grid_wide   = linspace(0.200, 0.400, 1000);                                 

% Pre-allocate storage (exact and approximate X for each grid / each regime)
[X_ex_nk, X_ap_nk] = deal(nan(size(kappa_grid_narrow)));  % narrow, vary κ
[X_ex_ns, X_ap_ns] = deal(nan(size(sigma_grid_narrow)));  % narrow, vary σ
[X_ex_wk, X_ap_wk] = deal(nan(size(kappa_grid_wide)));    % wide,   vary κ
[X_ex_ws, X_ap_ws] = deal(nan(size(sigma_grid_wide)));    % wide,   vary σ

% Exact first-order condition (eq. 16) as a function of X for fixed (κ,σ)
foc_exact = @(X, kappa, sigma) ...
    (2*X) / (3*beta^2) ...
    - 5*tanh(sqrt(2*beta/sigma^2)*X) ./ (2*(sqrt(2*beta/sigma^2))^3 .* X.^2 .* beta^2) ...
    + 5*(sqrt(2*beta/sigma^2)) ./ (2*(sqrt(2*beta/sigma^2))^3 .* X .* beta^2 ...
        .* cosh(sqrt(2*beta/sigma^2)*X).^2) ...
    + sinh(2*sqrt(2*beta/sigma^2)*X) ./ (2*(sqrt(2*beta/sigma^2)) .* beta^2 ...
        .* cosh(sqrt(2*beta/sigma^2)*X).^4) ...
    - (kappa * sigma^2) ./ (2 * X.^2);

% Narrow-band approximated FOC: d/dX [ (17/315)(η^4 X^6/β^2) + κσ²/(2X) ] = 0
foc_narrow = @(X, kappa, sigma) ...
    (408/315) * X.^5 / sigma^4 - (kappa * sigma^2) / (2 * X.^2);

% Wide-band approximated FOC: d/dX [ X²/(3β²) + κσ²/(2X) ] = 0
foc_wide = @(X, kappa, sigma) ...
    2 * X / (3 * beta^2) - (kappa * sigma^2) / (2 * X.^2);

% Helper: b from X for wide-band regime (b ≈ X/β)
b_wide_approx = @(X) X / beta;

% Helper: b from X for narrow-band regime (b ≈ 2X³/(3σ²))
b_narrow_approx = @(X, sigma) 2 * X^3 / (3 * sigma^2);

% Helper: exact b from X and σ
b_exact = @(X, sigma) ...
    (1 / (sqrt(2*beta/sigma^2) * beta)) * ...
    (sqrt(2*beta/sigma^2) * X - tanh(sqrt(2*beta/sigma^2) * X));

% --- Varying κ (narrow-band regime) ---
for i = 1:length(kappa_grid_narrow)
    kappa_i = kappa_grid_narrow(i);
    sigma_i = sigma_fix_narrow;
    try
        X_ex_nk(i) = fzero(@(X) foc_exact(X, kappa_i, sigma_i),  [1e-4, 5]);
        X_ap_nk(i) = fzero(@(X) foc_narrow(X, kappa_i, sigma_i), [1e-4, 5]);
    catch; end
end

% --- Varying σ (narrow-band regime) ---
for i = 1:length(sigma_grid_narrow)
    kappa_i = kappa_fix_narrow;
    sigma_i = sigma_grid_narrow(i);
    try
        X_ex_ns(i) = fzero(@(X) foc_exact(X, kappa_i, sigma_i),  [1e-4, 5]);
        X_ap_ns(i) = fzero(@(X) foc_narrow(X, kappa_i, sigma_i), [1e-4, 5]);
    catch; end
end

% --- Varying κ (wide-band regime) ---
for i = 1:length(kappa_grid_wide)
    kappa_i = kappa_grid_wide(i);
    sigma_i = sigma_fix_wide;
    try
        X_ex_wk(i) = fzero(@(X) foc_exact(X, kappa_i, sigma_i), [1e-4, 5]);
        X_ap_wk(i) = fzero(@(X) foc_wide(X, kappa_i, sigma_i),  [1e-4, 5]);
    catch; end
end

% --- Varying σ (wide-band regime) ---
for i = 1:length(sigma_grid_wide)
    kappa_i = kappa_fix_wide;
    sigma_i = sigma_grid_wide(i);
    try
        X_ex_ws(i) = fzero(@(X) foc_exact(X, kappa_i, sigma_i), [1e-4, 5]);
        X_ap_ws(i) = fzero(@(X) foc_wide(X, kappa_i, sigma_i),  [1e-4, 5]);
    catch; end
end

% Log deviations (in percent)
ld_nk = 100 * log(X_ap_nk ./ X_ex_nk);
ld_ns = 100 * log(X_ap_ns ./ X_ex_ns);
ld_wk = 100 * log(X_ap_wk ./ X_ex_wk);
ld_ws = 100 * log(X_ap_ws ./ X_ex_ws);

fig8 = figure('Position', [100 100 1400 700], 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

ax1 = nexttile;
hold on;
plot(kappa_grid_narrow, ld_nk, 'LineWidth', 2.5, 'Color', clr_approx_small);
yline(0, '--k', 'LineWidth', 1.5);
xline(kappa_fix_narrow, ':', 'Color', clr_baseline, 'LineWidth', 2.0, ...
    'Label', 'Baseline', 'Interpreter', 'latex', 'FontSize', 16);
hold off;
xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex');
title(sprintf('Varying $\\kappa$  ($\\sigma = %.3f$)', sigma_fix_narrow), ...
    'Interpreter', 'latex', 'FontSize', 20);
xlim([min(kappa_grid_narrow), max(kappa_grid_narrow)]);
ylim([-10, 10]); grid on; grid minor;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 16);

ax2 = nexttile;
hold on;
plot(sigma_grid_narrow, ld_ns, 'LineWidth', 2.5, 'Color', clr_approx_small);
yline(0, '--k', 'LineWidth', 1.5);
xline(sigma_fix_narrow, ':', 'Color', clr_baseline, 'LineWidth', 2.0, ...
    'Label', 'Baseline', 'Interpreter', 'latex', 'FontSize', 16);
hold off;
xlabel('$\sigma$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex');
title(sprintf('Varying $\\sigma$  ($\\kappa = %.3f$)', kappa_fix_narrow), ...
    'Interpreter', 'latex', 'FontSize', 20);
xlim([min(sigma_grid_narrow), max(sigma_grid_narrow)]);
ylim([-10, 10]); grid on; grid minor;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 16);

ax3 = nexttile;
hold on;
plot(kappa_grid_wide, ld_wk, 'LineWidth', 2.5, 'Color', clr_approx_large);
yline(0, '--k', 'LineWidth', 1.5);
xline(kappa_fix_wide, ':', 'Color', clr_baseline, 'LineWidth', 2.0, ...
    'Label', 'Baseline', 'Interpreter', 'latex', 'FontSize', 16);
hold off;
xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex');
title(sprintf('Varying $\\kappa$  ($\\sigma = %.3f$)', sigma_fix_wide), ...
    'Interpreter', 'latex', 'FontSize', 20);
xlim([min(kappa_grid_wide), max(kappa_grid_wide)]);
ylim([-10, 10]); grid on; grid minor;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 16);

ax4 = nexttile;
hold on;
plot(sigma_grid_wide, ld_ws, 'LineWidth', 2.5, 'Color', clr_approx_large);
yline(0, '--k', 'LineWidth', 1.5);
xline(sigma_fix_wide, ':', 'Color', clr_baseline, 'LineWidth', 2.0, ...
    'Label', 'Baseline', 'Interpreter', 'latex', 'FontSize', 16);
hold off;
xlabel('$\sigma$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Log deviation (approx / exact) $\times 100$, \%', 'Interpreter', 'latex');
title(sprintf('Varying $\\sigma$  ($\\kappa = %.3f$)', kappa_fix_wide), ...
    'Interpreter', 'latex', 'FontSize', 20);
xlim tight; ylim([-10, 10]); grid on; grid minor;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridAlpha', 0.15, 'FontSize', 16);

% Row labels (left-side annotations)
ax_lbl = axes(fig8, 'Position', [0 0 1 1], 'Visible', 'off');
text(ax_lbl, 0.015, 0.75, 'Narrow band approximation', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');
text(ax_lbl, 0.015, 0.25, 'Wide band approximation', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');

% =========================================================================
%  SECTION 3 — AUXILIARY MODEL FUNCTIONS  (Appendix Figures 8 - 10)
%
%  Plots q(v), dq(v), h(v), h²(v), dh(v) — the functions used in the
%  calibration algorithm.  These confirm the monotonicity properties
%  that guarantee a unique solution to the inverse calibration problem.
% =========================================================================

fprintf('\n=== SECTION 3 : Auxiliary functions q(v) and h(v) ===\n');

v_grid = linspace(1e-4, 5, 1000);

h_vals  = h_func(v_grid, beta);
dh_vals = dh_func(v_grid, beta);
q_vals  = q_func(v_grid, beta);
dq_vals = dq_func(v_grid, beta);

% Appendix Figure 7: h(v), h²(v), dh(v)
figure('Color', 'w');
subplot(3, 1, 1);
plot(v_grid, h_vals,    'LineWidth', 2, 'Color', COL.band15);
title('$h(v)$',  'Interpreter', 'latex', 'FontSize', 20); grid on;

subplot(3, 1, 2);
plot(v_grid, h_vals.^2, 'LineWidth', 2, 'Color', COL.band15);
title('$h(v)^2$', 'Interpreter', 'latex', 'FontSize', 20); grid on;

subplot(3, 1, 3);
plot(v_grid, dh_vals,   'LineWidth', 2, 'Color', COL.band15);
title('$dh(v)$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$v$',    'Interpreter', 'latex', 'FontSize', 20); grid on;

% Appendix Figure 6: q(v), dq(v)
figure('Color', 'w');
subplot(2, 1, 1);
plot(v_grid, q_vals,  'LineWidth', 2, 'Color', COL.band15);
title('$q(v)$',  'Interpreter', 'latex', 'FontSize', 20); grid on;

subplot(2, 1, 2);
plot(v_grid, dq_vals, 'LineWidth', 2, 'Color', COL.band15);
title('$dq(v)$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$v$',    'Interpreter', 'latex', 'FontSize', 20); grid on;

% -------------------------------------------------------------------------
%  psi(t): second-order-condition check for the exact FOC (footnote in
%  "Closed-Form Solutions" / Proposition 3)
% -------------------------------------------------------------------------
fprintf('\n=== Verifying psi(t) > 0 for t > 0 (exact-case S.O.C.) ===\n');

t_grid = linspace(1e-6, 15, 5000);   % t=0 excluded: tanh(0)/0 is 0/0
psi_vals = psi_func(t_grid);

fprintf('min psi(t) on grid: %.10f (at t = %.6f)\n', min(psi_vals), t_grid(psi_vals == min(psi_vals)));
fprintf('psi(t) as t->0+:  %.8f, %.8f, %.8f  (t = 1e-6, 1e-4, 1e-2)\n', ...
    psi_func(1e-6), psi_func(1e-4), psi_func(1e-2));
fprintf('psi(t) as t->inf: %.8f, %.8f  (t = 10, 50)\n', psi_func(10), psi_func(50));
fprintf('Any non-positive values found? %s\n', mat2str(any(psi_vals <= 0)));

figure('Color', 'w');
plot(t_grid, psi_vals, 'LineWidth', 2, 'Color', COL.band15);
yline(0, '--k', 'LineWidth', 1);
title('$\psi(t)$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\psi(t)$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$',    'Interpreter', 'latex', 'FontSize', 16); 
xlim([0, 15]); ylim([-0.2, 2.2]); grid on;


%% =========================================================================
%  SECTION 4 — ITALIAN LIRA / ECU DATA  (Figures 3 & 4)
%
%  We read two ECB CSV files:
%   - The main file covers the bulk of the ERM I history.
%   - The supplementary file covers the tail of the series (appended below).
%  Both are concatenated and trimmed to the ERM I window 1979–1998.
%
%  The estimation window is 12 Jan 1987 – 16 Sep 1992:
%    ±6 % sub-period : 12 Jan 1987 – 4  Jan 1990
%    ±2.25% sub-period:  5 Jan 1990 – 11 Sep 1992
% =========================================================================

fprintf('\n=== SECTION 4 : Loading Italian lira / ECU data ===\n');

% Read and concatenate the two ECB data files
% The files are ECB Statistical Data Warehouse exports; columns G:H
% contain (TIME_PERIOD, OBS_VALUE).
first_part  = readtable(fullfile(DATA_DIR, 'ert_h_eur_d__custom_19839545_linear.csv'), ...
    'Range', 'G:H');
second_part = readtable(fullfile(DATA_DIR, 'ert_h_eur_d__custom_19839614_linear.csv'), ...
    'Range', 'G1:H288');

data_raw = [first_part; second_part];
data_raw = data_raw(240:end, :);   % trim to start at the ERM entry (Mar 1979)

t = data_raw.TIME_PERIOD;   % datetime vector (weekdays only)
s = data_raw.OBS_VALUE;     % ITL per ECU (nominal exchange rate level)

% -------------------------------------------------------------------------
%  ERM I central parities and realignment dates (Table 2 in thesis)
%  Source: Official EU/Commission records (EC Bulletins).
% -------------------------------------------------------------------------
erm6 = table( ...
    [datetime(1979,3,12); datetime(1979,9,24); datetime(1981,3,23); ...
     datetime(1981,10,5); datetime(1982,6,14); datetime(1983,3,21); ...
     datetime(1985,7,22); datetime(1986,4,7);  datetime(1987,1,12)], ...
    [datetime(1979,9,24); datetime(1981,3,23); datetime(1981,10,5); ...
     datetime(1982,6,14); datetime(1983,3,21); datetime(1985,7,22); ...
     datetime(1986,4,7);  datetime(1987,1,12); datetime(1990,1,5)], ...
    [1148.15; 1159.42; 1262.92; 1300.67; 1350.27; ...
     1386.78; 1520.06; 1496.21; 1483.58], ...
    'VariableNames', {'Start', 'End', 'Parity'});

% -------------------------------------------------------------------------
%  Figure 3 — Full ERM I history (1979–1998)
% -------------------------------------------------------------------------
fig3 = figure('Position', [100 100 1600 900], 'Color', 'w');
hold on; grid on;
plot(t, s, '-', 'LineWidth', 1.7, 'Color', COL.obs, 'DisplayName', 'Observed ITL/ECU');

% Draw each ±6 % sub-period as a shaded band with central-parity dotted line
band6_pct = 0.06;
for i = 1:height(erm6)
    ref = erm6.Parity(i);
    fill([erm6.Start(i) erm6.End(i) erm6.End(i) erm6.Start(i)], ...
         ref * [1-band6_pct 1-band6_pct 1+band6_pct 1+band6_pct], ...
         COL.band6, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot([erm6.Start(i) erm6.End(i)], [ref ref], ':', ...
         'Color', COL.band6, 'LineWidth', 2.0, 'HandleVisibility', 'off');
    xline(erm6.Start(i), 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
end
% Single legend entry for all ±6 % periods
patch(nan, nan, COL.band6, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
    'DisplayName', 'ERM ±6% band');

% ±2.25 % band  (Jan 1990 – Sep 1992)
ref2   = 1529.70;  band2_pct = 0.0225;
p2_start = datetime(1990,1,5);  p2_end = datetime(1992,9,16);
fill([p2_start p2_end p2_end p2_start], ...
     ref2 * [1-band2_pct 1-band2_pct 1+band2_pct 1+band2_pct], ...
     COL.band2, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'ERM ±2.25% band');
plot([p2_start p2_end], [ref2 ref2], ':', 'Color', COL.band2, ...
     'LineWidth', 2.0, 'HandleVisibility', 'off');

% ERM suspension (Sep 1992 – Nov 1996)
ns = datetime(1992,9,16);  ne = datetime(1996,11,24);
yl = ylim;
fill([ns ne ne ns], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.85 0.85 0.85], 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', 'Outside ERM');

% ±15 % band  (Nov 1996 – Jan 1999)
ref3   = 1906.48;  band15_pct = 0.15;
p3_start = datetime(1996,11,24);  p3_end = datetime(1999,1,1);
fill([p3_start p3_end p3_end p3_start], ...
     ref3 * [1-band15_pct 1-band15_pct 1+band15_pct 1+band15_pct], ...
     COL.band15, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'ERM ±15% band');
plot([p3_start p3_end], [ref3 ref3], ':', 'Color', COL.band15, ...
     'LineWidth', 2.0, 'HandleVisibility', 'off');

xlabel('Time', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Italian lira per ECU', 'FontSize', 14, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 18;  ax.Box = 'on';  ax.GridAlpha = 0.25;
ax.XAxis.TickLabelFormat = 'MMM yyyy';
ax.XAxis.TickValues = linspace(min(t), max(t), 18);
xtickangle(40);
legend('Location', 'northwest', 'FontSize', 18);
uistack(findobj(gca, 'DisplayName', 'Observed ITL/ECU'), 'top');
hold off;

% -------------------------------------------------------------------------
%  Figure 4 — Estimation window: Jan 1987 – Sep 1992
% -------------------------------------------------------------------------

% Estimation window endpoints
t_est_start = datetime(1987, 1, 12);
t_est_end   = datetime(1992, 9, 16);

% Sub-period boundaries (also used for calibration below)
p6_start = datetime(1987, 1, 12);  p6_end = datetime(1990, 1, 4);   ref6 = 1483.58;
p2_start = datetime(1990, 1, 5);   p2_end = datetime(1992, 9, 16);  ref2 = 1529.70;

idx_est = t >= t_est_start & t <= t_est_end;

fig4 = figure('Position', [100 100 1600 900], 'Color', 'w');
hold on; grid on;

plot(t(idx_est), s(idx_est), '-', 'Color', COL.obs, 'LineWidth', 1.7, ...
    'DisplayName', 'Observed ITL/ECU');

% ±6 % band
fill([p6_start p6_end p6_end p6_start], ...
     ref6 * [1-band6_pct 1-band6_pct 1+band6_pct 1+band6_pct], ...
     COL.band6, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'ERM ±6% band');
plot([p6_start p6_end], [ref6 ref6], ':', 'Color', COL.band6, ...
     'LineWidth', 2.2, 'HandleVisibility', 'off');
xline(p6_start, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');

% ±2.25 % band
fill([p2_start p2_end p2_end p2_start], ...
     ref2 * [1-band2_pct 1-band2_pct 1+band2_pct 1+band2_pct], ...
     COL.band2, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'ERM ±2.25% band');
plot([p2_start p2_end], [ref2 ref2], ':', 'Color', COL.band2, ...
     'LineWidth', 2.2, 'HandleVisibility', 'off');
xline(p2_start, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');

xlim([t_est_start t_est_end]);
xlabel('Time', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Italian lira per ECU', 'FontSize', 13, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 18;  ax.GridAlpha = 0.25;  ax.Box = 'on';
ax.XAxis.TickLabelFormat = 'MMM yyyy';
ax.XAxis.TickValues = linspace(t_est_start, t_est_end, 12);
xtickangle(40);
legend('Location', 'northwest', 'FontSize', 18);
uistack(findobj(gca, 'DisplayName', 'Observed ITL/ECU'), 'top');
hold off;


% =========================================================================
%  ITALIAN LIRA / ECU DATA DEVIATION FROM CENTRAL PARITY (Figures 3 & 4)
%
%  Figures 6 and 7 in Appendix plot the lira's DEVIATION FROM ITS PREVAILING CENTRAL
%  PARITY (in percent), rather than the raw exchange-rate level. This is
%  deliberate: every claim made elsewhere in the thesis is about the lira's
%  position within its band, and a levels plot cannot show "pressure
%  against the depreciation boundary", because the boundary itself (the
%  band around the central parity) moves at every realignment. Plotting
%  the deviation makes the within-band position -- and the credibility
%  issue of frequent parity resets -- directly visible in the same
%  picture, with band limits drawn as shaded, color-filled step regions
%  that jump at each realignment, and every realignment date marked.
% =========================================================================

fprintf('\n=== SECTION 4 : Loading Italian lira / ECU data ===\n');

% Read and concatenate the two ECB data files
% The files are ECB Statistical Data Warehouse exports; columns G:H
% contain (TIME_PERIOD, OBS_VALUE).
first_part  = readtable(fullfile(DATA_DIR, 'ert_h_eur_d__custom_19839545_linear.csv'), ...
    'Range', 'G:H');
second_part = readtable(fullfile(DATA_DIR, 'ert_h_eur_d__custom_19839614_linear.csv'), ...
    'Range', 'G1:H288');

data_raw = [first_part; second_part];
data_raw = data_raw(240:end, :);   % trim to start at the ERM entry (Mar 1979)

t = data_raw.TIME_PERIOD;   % datetime vector (weekdays only)
s = data_raw.OBS_VALUE;     % ITL per ECU (nominal exchange rate level)

% -------------------------------------------------------------------------
%  ERM I central parities and realignment dates (Table 2 in thesis)
%  Source: Official EU/Commission records (EC Bulletins).
% -------------------------------------------------------------------------
erm6 = table( ...
    [datetime(1979,3,12); datetime(1979,9,24); datetime(1981,3,23); ...
     datetime(1981,10,5); datetime(1982,6,14); datetime(1983,3,21); ...
     datetime(1985,7,22); datetime(1986,4,7);  datetime(1987,1,12)], ...
    [datetime(1979,9,24); datetime(1981,3,23); datetime(1981,10,5); ...
     datetime(1982,6,14); datetime(1983,3,21); datetime(1985,7,22); ...
     datetime(1986,4,7);  datetime(1987,1,12); datetime(1990,1,5)], ...
    [1148.15; 1159.42; 1262.92; 1300.67; 1350.27; ...
     1386.78; 1520.06; 1496.21; 1483.58], ...
    'VariableNames', {'Start', 'End', 'Parity'});

% -------------------------------------------------------------------------
%  Sub-period boundaries used for calibration (moved up from Section 5,
%  since Figure 4 below needs them). Defined here once; Section 5 should
%  reuse these rather than redefining them.
% -------------------------------------------------------------------------
p6_start = datetime(1987, 1, 12);  p6_end = datetime(1990, 1, 4);   ref6 = 1483.58;
p2_start = datetime(1990, 1, 5);   p2_end = datetime(1992, 9, 16);  ref2 = 1529.70;

% -------------------------------------------------------------------------
%  Build a single "regimes" table covering the ENTIRE sample: each row is
%  a period with a constant prevailing parity and band width. This extends
%  erm6 (the ±6% sub-periods) with the ±2.25%, outside-ERM (no band), and
%  ±15% regimes.
% -------------------------------------------------------------------------
regimes = erm6(:, {'Start','End','Parity'});
regimes.BandPct = repmat(0.06, height(regimes), 1);
regimes.Label   = repmat({'±6%'}, height(regimes), 1);

regimes = [regimes; { ...
    datetime(1990,1,5), datetime(1992,9,16), 1529.70, 0.0225, {'±2.25%'} }];

regimes = [regimes; { ...
    datetime(1992,9,16), datetime(1996,11,24), NaN, NaN, {'Outside ERM'} }];

regimes = [regimes; { ...
    datetime(1996,11,24), max(t), 1906.48, 0.15, {'±15%'} }];

% -------------------------------------------------------------------------
%  Assign each observation to its prevailing parity and band width, and
%  compute the percentage deviation from that parity. NaN during the
%  outside-ERM float (no central parity, no band -> plotted as a gap).
% -------------------------------------------------------------------------
parity_t = nan(size(t));
band_t   = nan(size(t));
for i = 1:height(regimes)
    idx = t >= regimes.Start(i) & t < regimes.End(i);
    parity_t(idx) = regimes.Parity(i);
    band_t(idx)   = regimes.BandPct(i);
end
idx_last = t == regimes.End(end);
parity_t(idx_last) = regimes.Parity(end);
band_t(idx_last)   = regimes.BandPct(end);

dev_pct = (s - parity_t) ./ parity_t * 100;   % deviation from parity, in %

% =========================================================================
%  Appendix Figure 6 — Full ERM I history, deviation from prevailing parity
% =========================================================================
fig3 = figure('Position', [100 100 1600 900], 'Color', 'w');
hold on; grid on;

% Shade the outside-ERM (floating) window first, so it sits behind the data
idx_outside = strcmp(regimes.Label, 'Outside ERM');
os = regimes.Start(idx_outside);  oe = regimes.End(idx_outside);
fill([os oe oe os], [-20 -20 20 20], [0.85 0.85 0.85], ...
     'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', 'Outside ERM');

% Band limits: filled region + dotted step-line edges, per regime row
% (skipped for Outside ERM, which has no defined band)
for i = 1:height(regimes)
    switch regimes.Label{i}
        case '±6%',    col = COL.band6;
        case '±2.25%', col = COL.band2;
        case '±15%',   col = COL.band15;
        otherwise, col = [];
    end
    if ~isempty(col)
        b = 100*regimes.BandPct(i);
        st = regimes.Start(i); en = regimes.End(i);
        fill([st en en st], [-b -b b b], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot([st en],  [b b], ':', 'Color', col, 'LineWidth', 2.0, 'HandleVisibility', 'off');
        plot([st en], -[b b], ':', 'Color', col, 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
end
% One legend entry per band width, shown as a filled swatch. Uses NaT
% (not nan) for the x-coordinate: this axes' x-ruler is already locked to
% datetime by the fill() calls above, and a plain numeric nan here would
% reset the ruler to numeric and break every datetime-based object already
% drawn (observed in testing -- do not replace NaT with nan on this axes).
fill(NaT(1,4), nan(1,4), COL.band6,  'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'ERM ±6% band');
fill(NaT(1,4), nan(1,4), COL.band2,  'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'ERM ±2.25% band');
fill(NaT(1,4), nan(1,4), COL.band15, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'ERM ±15% band');

% Zero line (the central parity itself)
yline(0, 'Color', COL.parity, 'LineWidth', 1.2, 'LineStyle', '-', ...
      'HandleVisibility', 'off');

% Observed deviation from the prevailing parity (gaps automatically where
% parity_t is NaN, i.e. during the outside-ERM float)
plot(t, dev_pct, '-', 'LineWidth', 1.3, 'Color', COL.obs, ...
     'DisplayName', 'Deviation from prevailing parity');

% Mark every realignment date (every regime-row Start except the very
% first, which is EMS entry rather than a realignment)
for i = 2:height(regimes)
    xline(regimes.Start(i), 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
end

xlabel('Time', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Deviation from prevailing central parity, %', 'FontSize', 14, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 18;  ax.Box = 'on';  ax.GridAlpha = 0.25;
ax.XAxis.TickLabelFormat = 'MMM yyyy';
ax.XAxis.TickValues = linspace(min(t), max(t), 18);
xtickangle(40);
ylim([-20, 20]);
legend('Location', 'northwest', 'FontSize', 16);
uistack(findobj(gca, 'DisplayName', 'Deviation from prevailing parity'), 'top');
hold off;

% =========================================================================
%  Appendix Figure 7 — Estimation window: Jan 1987 – Sep 1992, deviation from
%  prevailing parity
% =========================================================================
t_est_start = datetime(1987, 1, 12);
t_est_end   = datetime(1992, 9, 16);
idx_est     = t >= t_est_start & t <= t_est_end;

dev_pct_est = dev_pct(idx_est);

fig4 = figure('Position', [100 100 1600 900], 'Color', 'w');
hold on; grid on;

% ±6% band: filled region + dotted step-line edges
fill([p6_start p6_end p6_end p6_start], [-6 -6 6 6], COL.band6, ...
     'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot([p6_start p6_end],  [6 6], ':', 'Color', COL.band6, 'LineWidth', 2.2, 'HandleVisibility', 'off');
plot([p6_start p6_end], -[6 6], ':', 'Color', COL.band6, 'LineWidth', 2.2, 'HandleVisibility', 'off');

% ±2.25% band: filled region + dotted step-line edges
fill([p2_start p2_end p2_end p2_start], [-2.25 -2.25 2.25 2.25], COL.band2, ...
     'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot([p2_start p2_end],  [2.25 2.25], ':', 'Color', COL.band2, 'LineWidth', 2.2, 'HandleVisibility', 'off');
plot([p2_start p2_end], -[2.25 2.25], ':', 'Color', COL.band2, 'LineWidth', 2.2, 'HandleVisibility', 'off');

% Legend swatches (NaT x-coordinate -- see note above Figure 3's equivalent)
fill(NaT(1,4), nan(1,4), COL.band6, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'ERM ±6% band');
fill(NaT(1,4), nan(1,4), COL.band2, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'ERM ±2.25% band');

yline(0, 'Color', COL.parity, 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(t(idx_est), dev_pct_est, '-', 'Color', COL.obs, 'LineWidth', 1.3, ...
     'DisplayName', 'Deviation from prevailing parity');

% Mark the single realignment inside this window (the Jan 1990 band change)
xline(p2_start, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

xlim([t_est_start t_est_end]);
xlabel('Time', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Deviation from prevailing central parity, %', 'FontSize', 13, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 18;  ax.GridAlpha = 0.25;  ax.Box = 'on';
ax.XAxis.TickLabelFormat = 'MMM yyyy';
ax.XAxis.TickValues = linspace(t_est_start, t_est_end, 12);
xtickangle(40);
ylim([-8, 8]);
legend('Location', 'northwest', 'FontSize', 18);
uistack(findobj(gca, 'DisplayName', 'Deviation from prevailing parity'), 'top');
hold off;


%% =========================================================================
%  SECTION 5 — BASELINE CALIBRATION  (Table 3)
%
%  Algorithm (Appendix A.2):
%   Step 1. Solve the ratio equation Var(s)/b² = q(v)/h(v)² for v = X/σ.
%           This eliminates σ and yields a unique v̂ (Newton in log space).
%   Step 2. Recover σ̂ = b_obs / h(v̂).
%   Step 3. Recover θ̂ = κ̂/σ̂ = g(v̂) from the first-order condition.
%   Step 4. Recover κ̂ = θ̂ · σ̂.
%
%  The procedure is applied independently to each band regime.
%
%  Inputs:
%   b_obs_6   = 0.06    (±6 % band half-width)
%   b_obs_2   = 0.0225  (±2.25 % band half-width)
%   β = 1/0.87 ≈ 1.1494 (Coenen & Vega, 2001)
% =========================================================================
 
fprintf('\n=== SECTION 5 : Baseline calibration ===\n');
 
% Exclude the three crisis days (14–16 Sep 1992) from the narrow-band sample
p2_end_est = datetime(1992, 9, 11);
 
% Compute fractional deviations of the exchange rate from its central parity:
%   x_t = (S_t - S̄) / S̄
idx_6 = t >= p6_start & t <= p6_end;
idx_2 = t >= p2_start & t <= p2_end_est;
 
x_6 = (s(idx_6) - ref6) / ref6;
x_2 = (s(idx_2) - ref2) / ref2;
 
% Sample variances — in decimal-squared units.
var_s_6 = mean((x_6 - mean(x_6)).^2);
var_s_2 = mean((x_2 - mean(x_2)).^2);
 
fprintf('Sample variance (±6%%  period): %.6f  [%.4f%%²]\n', var_s_6, var_s_6*100^2);
fprintf('Sample variance (±2.25%% period): %.6f  [%.4f%%²]\n', var_s_2, var_s_2*100^2);
 
b_obs_6 = 0.06;      
b_obs_2 = 0.0225;   
 
% Run inverse calibration for each regime
fprintf('\n--- Calibrating ±6%% regime ---\n');
res_6pct = inverse_calibrate(var_s_6, b_obs_6, beta);
 
fprintf('\n--- Calibrating ±2.25%% regime ---\n');
res_2pct = inverse_calibrate(var_s_2, b_obs_2, beta);

% -------------------------------------------------------------------------
%  Print Table 3
% -------------------------------------------------------------------------
fprintf('\n=== TABLE 3 : Calibration results ===\n');
fprintf('%-30s %12s %12s\n', '', '±6% regime', '±2.25% regime');
fprintf('%-30s %12s %12s\n', '', '(1987–1990)', '(1990–1992)');
fprintf('%s\n', repmat('-', 1, 56));
fprintf('%-30s %12.2f%% %11.2f%%\n',  'Observed band b_obs',      b_obs_6*100,           b_obs_2*100);
fprintf('%-30s %11.4f%%%% %10.4f%%%%\n', 'Observed variance Var(s)', var_s_6*100^2,  var_s_2*100^2);
fprintf('%-30s %12.4f %12.4f\n',      'Estimated X*',              res_6pct.X*100,        res_2pct.X*100);
fprintf('%-30s %12.4f %12.4f\n',      'Estimated sigma',           res_6pct.sigma,    res_2pct.sigma);
fprintf('%-30s %12.4f %12.4f\n',      'Estimated kappa',           res_6pct.kappa,    res_2pct.kappa);
fprintf('%-30s %12.4f %12.4f\n',      'v = X*/sigma',              res_6pct.v,        res_2pct.v);
fprintf('%-30s %12.4f %12.4f\n',      'theta = kappa/sigma',       res_6pct.theta,    res_2pct.theta);
 
%% =========================================================================
%  SECTION 6 — COUNTERFACTUAL ANALYSIS  (Table 4)
%
%  We hold one structural parameter fixed at its estimated value in one
%  regime and substitute the other regime's estimate, then re-solve the
%  planner's optimality condition to get the counterfactual band.
%
%  CF1: (κ̂_6%, σ̂_2.25%)  – isolates the volatility channel
%  CF2: (κ̂_2.25%, σ̂_6%)  – isolates the cost channel
% =========================================================================
 
fprintf('\n=== SECTION 6 : Counterfactual exercise ===\n');
 
% CF1: ±6% intervention cost + ±2.25% fundamental volatility
cf1 = forward_solve_optimal_band(res_6pct.kappa, res_2pct.sigma, beta);
fprintf('CF1 — kappa_6%%, sigma_2.25%%:\n');
fprintf('   theta = %.4f,  v* = %.4f,  X* = %.4f,  b* = %.4f\n', ...
    res_6pct.kappa / res_2pct.sigma, cf1.v, cf1.X, cf1.b);
 
% CF2: ±2.25% intervention cost + ±6% fundamental volatility
cf2 = forward_solve_optimal_band(res_2pct.kappa, res_6pct.sigma, beta);
fprintf('CF2 — kappa_2.25%%, sigma_6%%:\n');
fprintf('   theta = %.4f,  v* = %.4f,  X* = %.4f,  b* = %.4f\n', ...
    res_2pct.kappa / res_6pct.sigma, cf2.v, cf2.X, cf2.b);
 
% -------------------------------------------------------------------------
%  Print Table 4
% -------------------------------------------------------------------------
fprintf('\n=== TABLE 4 : Baseline and counterfactual optimal bands ===\n');
fprintf('%-28s %10s %10s %10s %10s\n', '', 'b=6%', 'b=2.25%', 'CF(k6,s2)', 'CF(k2,s6)');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%-28s %10.2f %10.2f %10s %10s\n', 'Observed band', ...
    b_obs_6*100, b_obs_2*100, '-', '-');
fprintf('%-28s %10.4f %10.4f %10s %10s\n', 'Observed Var(s) (%^2)', ...
    var_s_6*100^2, var_s_2*100^2, '-', '-');
fprintf('%-28s %10.4f %10.4f %10.4f %10.4f\n', 'sigma', ...
    res_6pct.sigma, res_2pct.sigma, res_2pct.sigma, res_6pct.sigma);
fprintf('%-28s %10.4f %10.4f %10.4f %10.4f\n', 'kappa', ...
    res_6pct.kappa, res_2pct.kappa, res_6pct.kappa, res_2pct.kappa);
fprintf('%-28s %10.4f %10.4f %10.4f %10.4f\n', 'theta = kappa/sigma', ...
    res_6pct.theta, res_2pct.theta, ...
    res_6pct.kappa/res_2pct.sigma, res_2pct.kappa/res_6pct.sigma);
fprintf('%-28s %10.4f %10.4f %10.4f %10.4f\n', 'v = X*/sigma', ...
    res_6pct.v, res_2pct.v, cf1.v, cf2.v);
fprintf('%-28s %10.4f %10.4f %10.4f %10.4f\n', 'X*', ...
    res_6pct.X*100, res_2pct.X*100, cf1.X*100, cf2.X*100);
 
% =========================================================================
%  SECTION 7 — FIGURE 5: Optimal X* vs σ for two κ values
%
%  Plots the optimal intervention trigger X*(σ) for the two calibrated
%  κ values, with vertical dashed lines at the estimated σ̂ values.
%  The four intersections correspond to the counterfactual entries in Table 4.
% =========================================================================
 
fprintf('\n=== SECTION 7 : Figure 5 — X*(σ) for two kappa values ===\n');
 
% -------------------------------------------------------------------------
%  Figure 5 — X*(σ) for two kappa values
% -------------------------------------------------------------------------
plot_Xstar_vs_sigma(beta, res_2pct.kappa, res_6pct.kappa, ...
    res_2pct.sigma, res_6pct.sigma, ...
    0.001, 0.15, 1000, COL.band15, COL.band6);

% =========================================================================
%  SECTION 8 — ROBUSTNESS TO β = 1/α  (Appendix Figures 9 & 10)
%
%  We re-run the calibration over a grid of money-demand semi-elasticities
%  α (equivalently β = 1/α) to assess sensitivity.
%  Baseline: α = 0.87, β ≈ 1.1494.
% =========================================================================
 
fprintf('\n=== SECTION 8 : Robustness over alpha (money-demand semi-elasticity) ===\n');
 
alpha_grid = [1/3, 0.50, 0.87, 1.00, 2.00];
beta_grid  = 1 ./ alpha_grid;
n_beta     = numel(beta_grid);
 
% Pre-allocate storage
fields = {'v', 'X', 'theta', 'sigma', 'kappa'};
rob_2 = struct(); rob_6 = struct();
for f = fields
    rob_2.(f{1}) = nan(n_beta, 1);
    rob_6.(f{1}) = nan(n_beta, 1);
end
 
fprintf('\n%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    'alpha', 'beta', 'regime', 'v', 'sigma', 'kappa', 'theta');
fprintf('%s\n', repmat('-', 1, 84));
 
for i = 1:n_beta
    beta_i  = beta_grid(i);
    alpha_i = alpha_grid(i);
 
    r6 = inverse_calibrate(var_s_6, b_obs_6, beta_i);
    r2 = inverse_calibrate(var_s_2, b_obs_2, beta_i);
 
    for f = fields
        rob_6.(f{1})(i) = r6.(f{1});
        rob_2.(f{1})(i) = r2.(f{1});
    end
 
    fprintf('%-12.2f %-12.4f %-12s %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
        alpha_i, beta_i, '±2.25%', r2.v, r2.sigma, r2.kappa, r2.theta);
    fprintf('%-12.2f %-12.4f %-12s %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
        alpha_i, beta_i, '±6%',    r6.v, r6.sigma, r6.kappa, r6.theta);
end
 
% -------------------------------------------------------------------------
%  Appendix Figure 13: calibrated κ and σ vs α
% -------------------------------------------------------------------------
figure('Name', 'Robustness: kappa and sigma vs alpha', 'Color', 'w');
 
subplot(1, 2, 1);
hold on;
plot(alpha_grid, rob_2.kappa, '-o', 'LineWidth', 3, 'Color', clr_approx_small, ...
    'DisplayName', 'Band 2.25%');
plot(alpha_grid, rob_6.kappa, '-s', 'LineWidth', 3, 'Color', clr_approx_large, ...
    'DisplayName', 'Band 6%');
xline(0.87, '--k', '$\alpha_{\rm base} = 0.87$', 'Interpreter', 'latex', ...
    'DisplayName', 'Baseline \alpha', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'LineWidth', 2, 'FontSize', 16);
xlabel('\alpha  (semi-elasticity)', 'FontSize', 16);
ylabel('\kappa', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 14);
title('Calibrated \kappa', 'FontSize', 16);
grid on; grid minor;
 
subplot(1, 2, 2);
hold on;
plot(alpha_grid, rob_2.sigma, '-o', 'LineWidth', 3, 'Color', clr_approx_small, ...
    'DisplayName', 'Band 2.25%');
plot(alpha_grid, rob_6.sigma, '-s', 'LineWidth', 3, 'Color', clr_approx_large, ...
    'DisplayName', 'Band 6%');
xline(0.87, '--k', '$\alpha_{\rm base} = 0.87$', 'Interpreter', 'latex', ...
    'DisplayName', 'Baseline \alpha', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'LineWidth', 2, 'FontSize', 16);
xlabel('\alpha  (semi-elasticity)', 'FontSize', 16);
ylabel('\sigma', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 14);
title('Calibrated \sigma', 'FontSize', 16);
grid on; grid minor;
 
% -------------------------------------------------------------------------
%  Appendix Figure 14: optimal X* vs β, both regimes
% -------------------------------------------------------------------------
figure('Name', 'Robustness: X* vs beta', 'Color', 'w');
hold on;
plot(beta_grid, rob_2.X, '-o', 'LineWidth', 3, 'Color', clr_approx_small, ...
    'DisplayName', '$X^\star_{2.25\%}$');
plot(beta_grid, rob_6.X, '-s', 'LineWidth', 3, 'Color', clr_approx_large, ...
    'DisplayName', '$X^\star_{6\%}$');
xline(1/0.87, '--k', '$\beta_{\rm base} = 1.1494$', ...
    'Interpreter', 'latex', 'FontSize', 16, 'LineWidth', 2, ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', ...
    'DisplayName', 'Baseline $\beta$');
xlabel('$\beta$  $(= 1/\alpha)$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$X^\star$',               'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 16);
grid on; grid minor;
ax = gca;
ax.MinorGridLineStyle = ':';  ax.GridAlpha = 0.25;  ax.MinorGridAlpha = 0.15;
ax.XMinorTick = 'on';  ax.YMinorTick = 'on';
set(gca, 'FontSize', 14);

 
%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

% -------------------------------------------------------------------------
%  psi(t)  —  second-order-condition function for the exact FOC (eq. 16)
%
%  psi(t) := 2 - sech^2(t) * [ 3*tanh(t)/t - sech^2(t) + 2*tanh^2(t) ]
%
%  From the proof of Proposition 3: the exact objective's second derivative,
%  after multiplying by X and substituting the F.O.C., reduces to
%  psi(t)/beta^2 with t := eta*X. psi(t) > 0 for all t > 0 is therefore the
%  second-order condition confirming the exact FOC's root is a minimum.
%  We do not pursue a fully analytic proof that psi(t) > 0 for all t > 0;
%  this function exists to verify the sign numerically instead.
% -------------------------------------------------------------------------
function val = psi_func(t)
    sech2_t = sech(t).^2;
    tanh_t  = tanh(t);
    val = 2 - sech2_t .* (3*tanh_t./t - sech2_t + 2*tanh_t.^2);
end


% -------------------------------------------------------------------------
%  h(v, β)  —  normalised boundary function  (eq. 20 in thesis)
%
%  The observed band b satisfies b = σ · h(v, β), where v = X/σ.
%  Derived from the boundary condition eq. (6): b = (1/ηβ)(ηX − tanh(ηX)).
%
%  h(v) = v/β − tanh(√(2β) v) / (β √(2β))
% -------------------------------------------------------------------------
function val = h_func(v, beta)
    sq2b = sqrt(2 * beta);
    val  = v ./ beta - tanh(sq2b .* v) ./ (beta .* sq2b);
end
 
% -------------------------------------------------------------------------
%  dh(v, β)  —  derivative of h with respect to v
%
%  dh(v)/dv = (1/β)(tanh(√(2β) v))²   =  (1/β) sech²(√(2β) v) ... simplified:
%  From h(v) = v/β − tanh(sq2b·v)/(β·sq2b):
%    dh/dv = 1/β − sech²(sq2b·v)/β  =  tanh²(sq2b·v)/β
% -------------------------------------------------------------------------
function val = dh_func(v, beta)
    sq2b = sqrt(2 * beta);
    val  = (1 / beta) .* tanh(sq2b .* v).^2;
end
 
% -------------------------------------------------------------------------
%  q(v, β)  —  normalised variance function  (eq. 21 in thesis)
%
%  The steady-state variance satisfies Var(s) = σ² · q(v, β).
%  Derived from eq. (10), after substituting X = σv, η = √(2β)/σ.
%
%  q(v) = v²/(3β²) − 1/β³ + 5 tanh(√(2β)v) / (4β³ √(2β) v)
%          − sech²(√(2β)v) / (4β³)
% -------------------------------------------------------------------------
function val = q_func(v, beta)
    sq2b = sqrt(2 * beta);
    val  = v.^2 ./ (3 * beta^2) ...
         - 1 / beta^3 ...
         + 5 .* tanh(sq2b .* v) ./ (4 * beta^3 .* sq2b .* v) ...
         - 1/(4 * beta^3) .* (1 - tanh(sq2b .* v).^2);
end
 
% -------------------------------------------------------------------------
%  dq(v, β)  —  derivative of q with respect to v
% -------------------------------------------------------------------------
function val = dq_func(v, beta)
    sq2b = sqrt(2 * beta);
    val  = 2 .* v ./ (3 * beta^2) ...
         + (5 ./ (4 * beta^3 .* sq2b)) .* ...
           (sq2b .* v .* sech(sq2b .* v).^2 - tanh(sq2b .* v)) ./ v.^2 ...
         + (sq2b ./ (2 * beta^3)) .* sech(sq2b .* v).^2 .* tanh(sq2b .* v);
end
 
% -------------------------------------------------------------------------
%  g(v, β)  —  normalised FOC function  (eq. 22 / Corollary 2)
%
%  At the optimum, g(v) = θ where θ = κ/σ.
%  g is strictly increasing, so the equation is uniquely invertible.
%
%  g(v) = 4v³/(3β²)
%         − 5 tanh(√(2β)v) / (2√2 β^{7/2})
%         + 5v sech²(√(2β)v) / (2β³)
%         + 2√2 v² tanh(√(2β)v) sech²(√(2β)v) / β^{5/2}
% -------------------------------------------------------------------------
function theta = g_func(v, beta)
    sq2b = sqrt(2 * beta);
    theta = (4 .* v.^3) ./ (3 .* beta^2) ...
          - 5 .* tanh(sq2b .* v) ./ (2 .* sqrt(2) .* beta^(7/2)) ...
          + (5.* v)./(2*beta^3*cosh(sq2b .* v).^2) ...
          + (sqrt(2)*v.^2 .* sinh(2*sqrt(2*beta)*v)) ./ (2*beta^(5/2)*cosh(sq2b.*v).^4);
end
 
% -------------------------------------------------------------------------
%  newton_method  —  Newton's method in log space  (Algorithm 2)
%
%  Solves f(v) = 0 using x = log(v) as the iteration variable.
%  Working in log space keeps v = exp(x) strictly positive throughout.
%
%  INPUTS
%    f       : objective as a function of v (not log v)
%    df      : derivative df/dv
%    x0      : initial guess in log space (x0 = log(v0))
%    tol     : convergence tolerance on |f(v)|
%    maxiter : maximum number of Newton steps
%
%  OUTPUT
%    root : converged value of v
% -------------------------------------------------------------------------
function root = newton_method(f, df, x0, tol, maxiter)
    x = x0;
    for iter = 1:maxiter
        v   = exp(x);
        fx  = f(v);
        dfx = df(v) * v;   % chain rule: d/dx f(e^x) = f'(e^x) · e^x
 
        % Safety check: if the Newton step in log-space is huge the iteration
        % is diverging, 1000 is just heuristics (robust)
        if abs(fx / dfx) > 1000
            error('newton_method: step too large (|f/df| > 1000). ');
        end
 
        x_new = x - fx / dfx;
 
        if abs(fx) < tol
            root = exp(x_new);
            return;
        end
        x = x_new;
    end
    error('newton_method: did not converge within %d iterations.', maxiter);
end
 
% -------------------------------------------------------------------------
%  inverse_calibrate  —  back out (σ̂, X̂, κ̂) from observables
%                         (Algorithm 1 in Appendix A.2)
%
%  INPUTS
%    var_data : empirical variance of the exchange rate within the regime
%    b_obs    : observed band half-width (as a decimal, e.g. 0.06)
%    beta     : money-demand parameter β = 1/α
%
%  ALGORITHM
%   1. Solve b_obs² q(v) − Var(s) h(v)² = 0 for v = X/σ using Newton
%      in log space (unique solution because q(v)/h(v)² is monotone).
%   2. Recover σ̂ = b_obs / h(v̂).
%   3. Recover θ̂ = g(v̂)  (from the FOC Corollary 2).
%   4. Recover κ̂ = θ̂ · σ̂.
%
%  OUTPUT (struct)
%    out.v     = v̂ = X̂/σ̂
%    out.X     = X̂ = optimal fundamental boundary
%    out.sigma = σ̂
%    out.kappa = κ̂
%    out.theta = θ̂ = κ̂/σ̂
% -------------------------------------------------------------------------
function out = inverse_calibrate(var_data, b_obs, beta)
 
    % Objective: b² q(v) − Var(s) h(v)² = 0
    obj  = @(v)  b_obs^2 * q_func(v, beta) - var_data * h_func(v, beta).^2;
    dobj = @(v)  b_obs^2 * dq_func(v, beta) - 2 * h_func(v, beta) .* var_data .* dh_func(v, beta);
 
    % Initial guess for v: from the wide-band approximation Var ≈ b²/3
    % v₀ ≈ b / √Var  (wide-band approximation; log-space Newton corrects from here)
    v0 = b_obs / sqrt(var_data);
    x0 = log(v0);
 
    % Solve in log space
    v_hat = exp(newton_method(obj, dobj, x0, 1e-20, 10000));
 
    % Recover structural parameters
    sigma_hat = sqrt(var_data / q_func(v_hat, beta));
    theta_hat = g_func(v_hat, beta);
    kappa_hat = theta_hat * sigma_hat;
    X_hat     = v_hat * sigma_hat;
 
    out.v     = v_hat;
    out.X     = X_hat;
    out.sigma = sigma_hat;
    out.theta = theta_hat;
    out.kappa = kappa_hat;
end
 
% -------------------------------------------------------------------------
%  newton_scalar  —  simple Newton's method in the original (v) space
%
%  INPUTS
%    f       : objective as a function of v
%    df      : derivative df/dv  (may be a finite-difference approximation)
%    v0      : initial guess
%    tol     : convergence tolerance on |f(v)|
%    maxiter : maximum Newton steps
%
%  OUTPUT
%    v : converged solution
% -------------------------------------------------------------------------
function v = newton_scalar(f, df, v0, tol, maxiter)
    v = v0;
    for iter = 1:maxiter
        fv  = f(v);
        dfv = df(v);
        v   = v - fv / dfv;
        if abs(f(v)) < tol
            return;
        end
    end
    error('newton_scalar: did not converge within %d iterations.', maxiter);
end
 
% -------------------------------------------------------------------------
%  solve_forward  —  solve g(v, β) = θ for v  (used inside forward_solve
%                    and plot_Xstar_vs_sigma)
%
%  Brackets the root by doubling the upper bound until a sign change is
%  found, then calls fzero.  This makes the solver robust to large θ values
%  where the root may lie far from v = 1.
%
%  INPUTS
%    theta_target : target value of θ = κ/σ
%    beta         : money-demand parameter β
%
%  OUTPUT
%    v_star : solution to g(v, β) = θ
% -------------------------------------------------------------------------
function v_star = solve_forward(theta_target, beta)
 
    obj  = @(v) g_func(v, beta) - theta_target;
 
    % Build bracket: expand upper bound until a sign change is found
    v_lo = 1e-6;
    v_hi = 10.0;
    for k = 1:50
        if obj(v_lo) * obj(v_hi) < 0
            break;
        end
        v_hi = v_hi * 2;
        if k == 50
            error('solve_forward: could not bracket root. theta_target = %.6f', ...
                  theta_target);
        end
    end
 
    v_star = fzero(obj, [v_lo, v_hi], optimset('TolX', 1e-12, 'Display', 'off'));
end
 
% -------------------------------------------------------------------------
%  forward_solve  —  compute optimal (v*, X*) for given (κ, σ)
%
%  Given a (κ, σ) pair, computes θ = κ/σ and solves the FOC g(v) = θ for
%  v*, then recovers X* = v* σ.  Used for the counterfactual exercise
%  (Section 6) and inside plot_Xstar_vs_sigma.
%
%  The g-function is defined inline here (matching the original code) to
%  make the function self-contained.  It is equivalent to g_func above.
%
%  INPUTS
%    kappa : intervention cost parameter κ
%    sigma : fundamental volatility σ
%    beta  : money-demand parameter β = 1/α
%
%  OUTPUT (struct)
%    out.v     = v* = X*/σ
%    out.X     = X* (optimal fundamental boundary)
%    out.theta = θ  = κ/σ
%    out.b     = b* (optimal band half-width, from b = σ h(v*))
% -------------------------------------------------------------------------

function out = forward_solve(kappa, sigma, beta)
    theta = kappa / sigma;
    
    % FOC: g(v, β) = θ solve for v using finite-difference Newton
    obj = @(v) g_func(v, beta) - theta;   
    
    % Finite-difference derivative of the objective (robust to expression complexity)
    dg = @(v) (obj(v * 1.0001) - obj(v * 0.9999)) / (v * 0.0002);
    
    % Initial guess v₀ = 1 (calibrated v ≈ 1 in both ERM I regimes)
    v0 = 1.0;
    v  = newton_scalar(obj, dg, v0, 1e-12, 1000);

    out.v     = v;
    out.theta = theta;
    out.X     = v * sigma;
    out.b     = h_func(v, beta) * sigma;
end
 
% -------------------------------------------------------------------------
%  forward_solve_optimal_band  —  thin wrapper kept for call-site clarity
%
%  Calls forward_solve and returns the same struct.  This name is used at
%  the two counterfactual call sites in Section 6, while solve_forward (the
%  scalar root-finder) is used inside plot_Xstar_vs_sigma.
% -------------------------------------------------------------------------
function out = forward_solve_optimal_band(kappa, sigma, beta)
    out = forward_solve(kappa, sigma, beta);
end
 
% -------------------------------------------------------------------------
%  plot_Xstar_vs_sigma  —  Figure 5 in thesis
%
%  Plots X*(σ) for two calibrated κ values over a range of σ.
%  Vertical dashed lines mark the two estimated σ̂ values.
%
%  INPUTS
%    beta         : money-demand parameter
%    kappa1       : first κ value  (here: κ̂_{2.25%})
%    kappa2       : second κ value (here: κ̂_{6%})
%    sigma1       : first estimated σ̂   (for dashed-line label)
%    sigma2       : second estimated σ̂  (for dashed-line label)
%    sigma_min    : lower end of σ grid
%    sigma_max    : upper end of σ grid
%    N            : number of grid points
%    color1, color2 : line colours
% -------------------------------------------------------------------------
function plot_Xstar_vs_sigma(beta, kappa1, kappa2, sigma1, sigma2, ...
                             sigma_min, sigma_max, N, color1, color2)
 
    sigma_grid = linspace(sigma_min, sigma_max, N);
    X1 = zeros(size(sigma_grid));
    X2 = zeros(size(sigma_grid));
 
    for i = 1:N
        sig = sigma_grid(i);
        v1 = solve_forward(kappa1 / sig, beta);   % solve g(v) = κ1/σ
        v2 = solve_forward(kappa2 / sig, beta);   % solve g(v) = κ2/σ
        X1(i) = v1 * sig;
        X2(i) = v2 * sig;
    end
 
    figure('Name', 'Optimal X vs sigma', 'NumberTitle', 'off', 'Color', 'w');
    plot(sigma_grid, X1, 'LineWidth', 2, 'Color', color1, ...
        'DisplayName', ['$\kappa_{2.25\%} = ', num2str(round(kappa1, 4), '%.4g'), '$']);
    hold on;
    plot(sigma_grid, X2, 'LineWidth', 2, 'Color', color2, ...
        'DisplayName', ['$\kappa_{6\%} = ',   num2str(round(kappa2, 4), '%.4g'), '$']);

    xline(sigma1, '--', ...
        ['$\hat{\sigma}_{2.25\%} = ', num2str(round(sigma1, 4), '%.4g'), '$'], ...
        'Interpreter', 'latex', 'FontSize', 18, ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(sigma2, '--', ...
        ['$\hat{\sigma}_{6\%} = ', num2str(round(sigma2, 4), '%.4g'), '$'], ...
        'Interpreter', 'latex', 'FontSize', 18, ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
 
    xlabel('$\sigma$',          'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$X^\star(\sigma)$', 'Interpreter', 'latex', 'FontSize', 20);
    legend('Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 16);
    grid on; grid minor;
    ax = gca;
    ax.MinorGridLineStyle = ':';  ax.GridAlpha = 0.25;  ax.MinorGridAlpha = 0.15;
    ax.XMinorTick = 'on';  ax.YMinorTick = 'on';
    set(gca, 'FontSize', 16);
end
 
