%% Cleaning
clear all, close all, clc, format short e;

%% Question 1
% Pen-and-paper proof is written in LaTeX report.

%% Question 2
% Approximate u'' for u(x)=cos(pi*x/2), n=100.

n = 100; 
h = 2 / (n+1);  % Spacing
x = -1 + (1:n)' * h; % internal vector point

% build matrix
% central dfference formula: u''(xi) ≈ (u(xi-1) - 2u(xi) + u(xi+1))/h^2
D = (diag(ones(n-1,1),1) + diag(-2*ones(n,1)) + diag(ones(n-1,1),-1)) * (1/h^2);

% build U
U = cos(pi*x/2);

% obtain the approximation
V = D*U;

% u(x) = cos(πx/2), u''(x) = -(π^2/4)cos(πx/2)
exact = -(pi^2/4)*cos(pi*x/2);

figure;
plot(x, exact, 'b-', x, V, 'ro');
legend('exact', 'approximation');
%% Question 3
% set n list that satisfies h half-decreasing, so according the O(h^2), the
% error should be reduced to 1/4
n_list = [4,9,19,39,79,159,319,639];
errors = zeros(size(n_list));
% expected_ratio = 1/4;

for i = 1:length(n_list)
    n = n_list(i); 
    h = 2 / (n+1);  % Spacing
    x = -1 + (1:n)' * h; % internal vector point
    
    % build matrix
    % central dfference formula: u''(xi) ≈ (u(xi-1) - 2u(xi) + u(xi+1))/h^2
    D = (diag(ones(n-1,1),1) + diag(-2*ones(n,1)) + diag(ones(n-1,1),-1)) * (1/h^2);
    
    % build U
    U = cos(pi*x/2);
    
    % obtain the approximation
    V = D*U;
    
    % u(x) = cos(πx/2), u''(x) = -(π^2/4)cos(πx/2)
    exact = -(pi^2/4)*cos(pi*x/2);

    % max error
    errors(i) = max(abs(V-exact));
end

h_list = 2./(n_list+1);
error_ratio = errors(1:end-1) ./ errors(2:end);
expected_ratio = (h_list(1:end-1)./h_list(2:end)).^2;

% plot
figure;
loglog(h_list, errors, 'bo-', 'LineWidth', 2);
hold on;
loglog(h_list, errors(1)* (h_list/h_list(1)).^2, 'r--', 'LineWidth', 2);
xlabel('h');
ylabel('Maximum Error');
legend('Numerical Error', 'O(h^2) Reference', 'Location', 'best');

% print
fprintf('h \t\t Error \t\t Error Ratio \t Excepted Ratio \n');
for i = 1:length(n_list)
    if i < length(n_list)
        fprintf('%.4f \t\t %.2e \t\t %.2f \t\t %.2f \n', h_list(i), errors(i), error_ratio(i),  expected_ratio(i));
    else
        fprintf('%.4f \t\t %.2e \t\t - \t\t - \n', h_list(i), errors(i));
    end
end

%% Question 4
n = 100; 
h = 2 / (n+1);  % Spacing
x = -1 + (1:n)' * h; % internal vector point

% build matrix
% central dfference formula: u''(xi) ≈ (u(xi-1) - 2u(xi) + u(xi+1))/h^2
D = (diag(ones(n-1,1),1) + diag(-2*ones(n,1)) + diag(ones(n-1,1),-1)) * (1/h^2);

% u''(x) = -f(x)
f = @(x) pi^2 * exp(cos(pi*x)) .* (cos(pi*x) - sin(pi*x).^2);
u_2 = -f(x);

% D * U = u_2, use backslach command
U = D \ u_2;

% u(x) = ecos(πx) − 1/e
excat = exp(cos(pi*x)) - 1/exp(1);

% error
error = max(abs(U - excat));
fprintf('Max error: %.2e\n', error);


% plot
figure;
plot(x, excat, 'b-', x, U, 'ro');
legend('Exact', 'Numerical');


%% Question 5
% Test the Forward Euler function
f = @(t, y) y;
tspan = [0, 1];
y0 = 1;
nsteps = 10;
[t, y] = ForwardEuler(f, tspan, y0, nsteps);

% exact result
exact = exp(t);

% plot
figure;
plot(t, y, 'bo-', t, exact, 'r-');
legend('Forward Eular', 'Exact');

%% Question 6
% set the ode
A = [0, 1; -5/2, -3];
f = @(t, y) A * y;

% inital condition
u0 = [1; -1/2];
t0 = 0;
T = 1;
tspan = [t0, T];

% set time step for test
nsteps_list = [5, 10, 20, 40, 80, 160, 320];
error = zeros(length(nsteps_list), 1);
% time step size
tss_list = zeros(length(nsteps_list), 1);

for i = 1:length(nsteps_list)
    nsteps = nsteps_list(i);
    tss = (T-t0) / nsteps;
    tss_list(i) = tss;

    % call ForwardEuler
    [tHist, uHist]= ForwardEuler(f, tspan, u0, nsteps);

    %calcu u(t) and v(t)
    u_exact = exp(-3*tHist/2) .* (cos(tHist/2) + 2* sin(tHist/2));
    v_exact = exp(-3*tHist/2) .* (-1/2*cos(tHist/2) - 7/2* sin(tHist/2));
    uHist_exact = [u_exact, v_exact];

    % calcu error
    error = max(max(abs(uHist - uHist_exact)));
    errors(i) = error;
end

% plot for u_exat
figure;
hold on;
plot(tHist, uHist(:, 1), 'b-', 'LineWidth', 1.5);
hold on
plot(tHist, u_exact, 'ro', 'MarkerSize', 3);
legend('Forward Eular', 'Exact');
title('u(t)')

% plot for v_exat
figure;
hold on;
plot(tHist, uHist(:, 2), 'b-', 'LineWidth', 1.5);
hold on
plot(tHist, v_exact, 'ro', 'MarkerSize', 3);
legend('Forward Eular', 'Exact');
title('v(t)')

%% Question 7
% Pen-and-paper proof is written in LaTeX report.

%% Question 8

%% Question 9

%% Question 10

%% Question 11

%% Question 12
