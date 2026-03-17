%% Cleaning
clear all, close all, clc, format short e;

%% Question 1
% See newton.m

%% Question 2
%  Use the Newton solver to compute the 10th root of 10, i.e., x^10 = 10,

% Define f and df for the equation f(x) = x^10 - 10 = 0
f  = @(x) x.^10 - 10;
df = @(x) 10 .* x.^9;

% Set initial guess and tolerance for Newton's method
% Because f(1) < 0 || f(1) > 0  => x between 1 and 0
x0 = 1.5;       % initial guess
tol = eps;      % tolerance for convergence
maxiter = 100;  % maximum number of iterations

% Call Newton
[x, res, xvec, resvec] = newton(f, df, x0, maxiter, tol);

% Verify using Matlab's built-in nthroot
xref = nthroot(10, 10);
abs_err = abs(x - xref);

% Report only what the assignment requires
fprintf('\nresult\n');
fprintf('  Newton approximation  = %.16f\n', x);
fprintf('  Reference nthroot     = %.16f\n', xref);
fprintf('  Absolute error        = %.3e\n', abs_err);
fprintf('  Machine precision ε_m = %.3e\n', eps);

% Check if error ≤ machine precision
if abs_err <= eps
    fprintf('  PASS: |x - nthroot| ≤ ε_m (requirement satisfied).\n\n');
else
    fprintf('  FAIL: |x - nthroot| > ε_m (not within machine accuracy).\n\n');
end

%% Question 3
% Demonstrate that newton exits with an error for inadmissible maxiter or tol
f  = @(x) x.^2 - 2;
df = @(x) 2*x;
x0 = 1.0;

%% Case A: maxiter not a positive integer
try
    maxiter = 0;     % invalid (must be positive integer)
    tol     = 1e-10; % valid
    newton(f, df, x0, maxiter, tol);
catch ME
    fprintf('Caught error (A): %s [%s]\n', ME.message, ME.identifier);
end

try
    maxiter = 10.5;  % invalid (must be integer)
    tol     = 1e-10; % valid
    newton(f, df, x0, maxiter, tol);
catch ME
    fprintf('Caught error (A2): %s [%s]\n', ME.message, ME.identifier);
end

%% Case B: tol not a positive finite scalar
try
    maxiter = 50;  % valid
    tol     = 0;   % invalid (must be > 0)
    newton(f, df, x0, maxiter, tol);
catch ME
    fprintf('Caught error (B): %s [%s]\n', ME.message, ME.identifier);
end

try
    maxiter = 50;      % valid
    tol     = NaN;     % invalid (must be finite)
    newton(f, df, x0, maxiter, tol);
catch ME
    fprintf('Caught error (B2): %s [%s]\n', ME.message, ME.identifier);
end

%% Question 4(a)
% Objective:
%   For γ = 3, consider the function
%       f(x) = cosh(x) + cos(x) - 3.
%   Show that Newton's method converges quadratically to its nonnegative root x*,
%   provided that the initial guess x0 is sufficiently close to x*.
%
% Step 1:
%   f(0) = cosh(0) + cos(0) - 3 = 1 + 1 - 3 = -1 < 0.
%   As x -> ∞,  cosh(x) -> ∞  =>  f(x) -> ∞.
%   Hence, by the intermediate value theorem, there exists at least one root x* > 0.
%   The second derivative f"(x) = cosh(x) - cos(x) satisfies:
%       f"(x) > 0  for all x > 0,
%   because  cosh(x) > 1  and  cos(x) ≤ 1  with equality only at x = 0.
%   Therefore, f is strictly convex on [0, ∞), implying that the zero is unique.
%
% Step 2:
%   f'(x) = sinh(x) - sin(x)
%   f"(x) = cosh(x) - cos(x)
%   Both derivatives exist and are continuous ⇒ f ∈ C^2(ℝ).
%   At the root x*, we have  f'(x*) = sinh(x*) - sin(x*) > 0.
%   Thus, f'(x*) ≠ 0.
%
% Step 3:
%   According to the theorem from the lecture notes:
%     If f ∈ C^2 in a neighborhood of x*,  f(x*) = 0,  and  f'(x*) ≠ 0,
%     then Newton's method
%         x_{k+1} = x_k - f(x_k)/f'(x_k)
%     converges quadratically to x*,  provided x0 is sufficiently close to x*.
%
% Conclusion:
%   All the above conditions are satisfied for
%       f(x) = cosh(x) + cos(x) - 3  (with γ = 3),
%   hence Newton's method converges quadratically to the unique nonnegative root x*.

%% Question 4(b)
% Numerical evidence of quadratic convergence for f(x)=cosh(x)+cos(x)-3 (γ=3)
% Using newton.m to show e_{k+1} ≈ C * e_k^2 asymptotically. 

% Define f and df
f  = @(x) cosh(x) + cos(x) - 3;
df = @(x) sinh(x) - sin(x);

% High-accuracy reference root x*
x0_ref   = 1.0;         % start near the nonnegative root
maxiterR = 100;
tolR     = eps^2;       % tighter than machine epsilon to make reference very accurate
[x_star, res_star] = newton(f, df, x0_ref, maxiterR, tolR);  % uses your Q1 function
fprintf('\n \reference root x* ≈ %.16f, |f(x*)|=%.3e\n', x_star, res_star);

% Run Newton again (normal tolerance) and collect error sequence e_k
x0      = 0.8;          % sufficiently close initial guess
maxiter = 20;
tol     = 1e-14;        % standard high accuracy
[~, ~, xvec, ~] = newton(f, df, x0, maxiter, tol);

% Error sequence e_k = |x_k - x*|
ek = abs(xvec - x_star);

% Empirical order p_k and constant C_k
% p_k = log( e_{k+1}/e_k ) / log( e_k / e_{k-1} )  (needs k>=2)
% C_k = e_{k+1} / (e_k^2)
pk = NaN(length(ek),1);
Ck = NaN(length(ek),1);

% Loop over iterations to compute p_k and C_k
for k = 2:length(ek)-1
    if ek(k)>0 && ek(k-1)>0
        pk(k) = log(ek(k+1)/ek(k)) / log(ek(k)/ek(k-1));
    end
    if ek(k)>0
        Ck(k) = ek(k+1) / (ek(k)^2);
    end
end

% Print a small convergence table (last few iterations)
fprintf('\nQ4(b): observed order and constants (last 6 usable rows)\n');
fprintf('   k     e_k         e_{k+1}     C_k = e_{k+1}/e_k^2    p_k\n');
fprintf(' -----------------------------------------------------------------------------\n');
kshow = max(2, length(ek)-6) : max(2, length(ek)-1);
for k = kshow
    fprintf(' %3d  % .3e    % .3e      % .3e         % .2f\n', ...
        k, ek(k), ek(k+1), Ck(k), pk(k));
end
fprintf('Note: Quadratic convergence <=> p_k -> 2 and C_k -> constant.\n');

% Optional: plot residuals or errors vs iteration (log scale)
figure;
semilogy(0:length(ek)-1, ek, 'o-','LineWidth',1); grid on;
xlabel('k'); ylabel('|x_k - x^*|');
title('Q4(b): Error decay (Newton)');

%% Question 4(c)
% Use newton to approximate x* for γ=2 with at most 20 iterations and tol=ε_m.

% Common settings
maxiter = 20;
tol     = eps;

% ----- γ = 2 -----
f2  = @(x) cosh(x) + cos(x) - 2;
df2 = @(x) sinh(x) - sin(x);
x0_2 = 0.1;                        % start near the nonnegative root
[~, ~, xvec2, resvec2] = newton(f2, df2, x0_2, maxiter, tol);

% ----- γ = 3 -----
f3  = @(x) cosh(x) + cos(x) - 3;
df3 = @(x) sinh(x) - sin(x);
x0_3 = 0.8;                        % as in Q4(b)
[~, ~, xvec3, resvec3] = newton(f3, df3, x0_3, maxiter, tol);

% Residual sequences {|f(x_k)|}
k2 = 0:numel(resvec2)-1;
k3 = 0:numel(resvec3)-1;

% Plot (log scale on y)
figure;
semilogy(k2, resvec2, 'o-','LineWidth',1); hold on;
semilogy(k3, resvec3, 's-','LineWidth',1); grid on;
xlabel('k'); ylabel('|f(x_k)|');
legend('\γ=2','\γ=3','Location','best');
title('Q4(c): Residuals vs iteration (Newton)');

% Report convergence status (per assignment requirement)
fprintf('γ=2: iterations used = %d, final |f(x)| = %.3e\n', numel(resvec2)-1, resvec2(end));
fprintf('γ=3: iterations used = %d, final |f(x)| = %.3e\n', numel(resvec3)-1, resvec3(end));

%%
% The plot shows the residuals |f(x_k)| versus iteration k for γ = 2 and γ = 3.
% For γ = 3 (orange curve), Initially the residuals first increase
% (because the initial values are slightly farther away) and then decrease
% rapidly, with the residuals dropping to ~1e-10 after 9 to 10
% iterations. It agrees with the theoretical analysis in Q4(a,b), because
% f'(x*) != 0 at the nonnegative root.
%
% For γ = 2 (blue curve), the residuals decrease much more slowly and
% approximately linearly on the logarithmic scale. This behavior confirms
% *linear convergence*, since the derivative at the root is zero (f'(0) =
% 0), violating the assumption required for quadratic convergence.
%
% Hence, the numerical results fully support the theoretical prediction:
% Newton's method converges quadratically when f'(x*) != 0 (γ = 3) and only
% linearly when f'(x*) = 0 (γ = 2).

%% Question 5
% See secant.m

%% Question 6
% Use secant to approximate a root of f(x) = x^10 - 10 with x0=4, x1=3, tol=1e-5,
% and provide numerical evidence that the method converges with order
%  q = (1 + √5)/2 ≈ 1.618034.

% Define function
f = @(x) x.^10 - 10;

% Parameters
x0 = 4; 
x1 = 3;
maxiter = 50;
tol = 1e-5;

% Call secant
[x_sec, res_sec, xvec, resvec] = secant(f, x0, x1, maxiter, tol);

% Reference root: x* = 10^(1/10)
xstar = nthroot(10, 10);

% Compute absolute errors e_k
ek = abs(xvec - xstar);

% Estimate the convergence order q using consecutive errors:
%   q_k = log(e_{k+1}/e_k) / log(e_k/e_{k-1})
n = length(ek);
fprintf('\n')
if n < 3
    fprintf('Not enough iterations to estimate the order (need at least 3 errors).\n');
else
    for k = 3 : n-1   % use n-1 to avoid accessing ek(k+1) out of range
        if ek(k-1) > 0 && ek(k) > 0 && ek(k+1) > 0
            qk = log(ek(k+1)/ek(k)) / log(ek(k)/ek(k-1));
            fprintf('k = %2d: q_k = %.6f\n', k, qk);
        end
    end
end

%% 
% The computed values of q_k converge to 1.618, confirming that the secant
% method converges with order q = (1 + √5)/2, as expected.