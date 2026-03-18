%% Cleaning
clear all, close all, clc, format short e;

%% Question 1
% Pen-and-paper proof is written in LaTeX report.

%% Question 2


%% Question 3


%% Question 4(a)


%% Question 4(b)
% Numerical evidence of quadratic convergence for f(x)=cosh(x)+cos(x)-3 (γ=3)
% Using newton.m to show e_{k+1} ≈ C * e_k^2 asymptotically. 


%% Question 4(c)

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

%% 
