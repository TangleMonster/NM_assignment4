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
spatial_errors = zeros(size(n_list));
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
    spatial_errors(i) = max(abs(V-exact));
end

h_list = 2./(n_list+1);
error_ratio = spatial_errors(1:end-1) ./ spatial_errors(2:end);
expected_ratio = (h_list(1:end-1)./h_list(2:end)).^2;

% plot
figure;
loglog(h_list, spatial_errors, 'bo-', 'LineWidth', 2);
hold on;
loglog(h_list, spatial_errors(1)* (h_list/h_list(1)).^2, 'r--', 'LineWidth', 2);
xlabel('h');
ylabel('Maximum Error');
legend('Numerical Error', 'O(h^2) Reference', 'Location', 'best');

% print
fprintf('h \t\t Error \t\t Error Ratio \t Excepted Ratio \n');
for i = 1:length(n_list)
    if i < length(n_list)
        fprintf('%.4f \t\t %.2e \t\t %.2f \t\t %.2f \n', h_list(i), spatial_errors(i), error_ratio(i),  expected_ratio(i));
    else
        fprintf('%.4f \t\t %.2e \t\t - \t\t - \n', h_list(i), spatial_errors(i));
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
% set the ode
A = [0, 1; -5/2, -3];
g = @(t) [0; 0];

% inital condition
u0 = [1; -1/2];
t0 = 0;
T = 1;
tspan = [t0, T];

% set time step for test
nsteps_list = [5, 10, 20, 40, 80, 160, 320];
% time step size
tss_list = zeros(length(nsteps_list), 1);

for i = 1:length(nsteps_list)
    nsteps = nsteps_list(i);
    tss = (T-t0) / nsteps;
    tss_list(i) = tss;

    % call ForwardEuler
    [tHist, uHist]= CrankNicholson(A, g, tspan, u0, nsteps);

    %calcu u(t) and v(t)
    u_exact = exp(-3*tHist/2) .* (cos(tHist/2) + 2* sin(tHist/2));
    v_exact = exp(-3*tHist/2) .* (-1/2*cos(tHist/2) - 7/2* sin(tHist/2));
    uHist_exact = [u_exact, v_exact];
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
%% Question 9
% Pen-and-paper proof is written in LaTeX report.

%% Question 10
% Pen-and-paper proof is written in LaTeX report.

%% Question 11
% spatial grid
n = 100; 
h = 1 / (n+1);  % Spacing
x = 0 + (1:n)' * h; % internal vector point (not including 0 and 1)
x_full = (0: h :1)';

% time grid
T=1;
nsteps = 100;
tss = T / nsteps;
t = (0:tss:T);

% build matrix
D = (diag(ones(n-1,1),1) + diag(-2*ones(n,1)) + diag(ones(n-1,1),-1)) * (1/h^2);

% Using Crank Nicholson
I = eye(n);
B = I - (tss/2)*D;
L = I + (tss/2)*D;

% initial condition
phi = sin(2*pi*x);
U = phi;
U_records = zeros(nsteps+1, n);
U_records(1, :) = phi';

for i = 1:nsteps
    U = B \ (L*U);
    U_records(i+1, :) = U';
end

% build the full matrix including the boundary points (0 and 1)
U_full = zeros(nsteps+1, n+2);
U_full(:, 1) = 0;  %left boundary
U_full(:, end) = 0;  % right boundary
U_full(:, 2:end-1) = U_records;

[X, T_mesh] = meshgrid(x_full, t);
U_exact = exp(-4*pi^2*T_mesh) .* sin(2*pi*X);

% plot
figure;
surf(X, T_mesh, U_full, 'EdgeColor', 'none');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Crank-Nicolson');

figure;
surf(X, T_mesh, U_exact, 'EdgeColor', 'none');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Analystical one');
%% Question 12
%% initial fixed time step settings
n_list = [5,10,20,40,80,160,320,640];
T = 0.05;
tss_fixed = 0.001;
nsteps_fixed = round(T / tss_fixed);

% pre_set
h_list = zeros(length(n_list), 1);
spatial_errors = zeros(length(n_list),1);

for i = 1:length(n_list)
    n = n_list(i);
    h = 1/(n+1);
    h_list(i) = h;

    x = (1:n)' * h;

    % build matrix
    D = (diag(ones(n-1,1),1) + diag(-2*ones(n,1)) + diag(ones(n-1,1),-1)) * (1/h^2);
    
    % Using Crank Nicholson
    I = eye(n);
    B = I - (tss_fixed/2)*D;
    L = I + (tss_fixed/2)*D;

    % initial condition
    phi = sin(2*pi*x);
    U = phi;

    % calcu the numerical solution for internal points at time t
    for step = 1:nsteps_fixed
        U = B \ (L*U);
    end

    % build the full matrix including the boundary points (0 and 1)
    x_full = (0:h:1)';
    U_full = zeros(1, n+2);
    U_full(1) = 0; %left boundary
    U_full(end) = 0; % right boundary
    U_full(2:end-1) = U';

    % excat soluction
    U_exact = exp(-4*pi^2*T) .* sin(2*pi*x_full);
    
    % calcu the error
    spatial_errors(i) = max(abs(U_full' - U_exact));
end 

% calcu the space convergence rates
space_error_rates = log(spatial_errors(1:end-1) ./ spatial_errors(2:end)) ./ log(h_list(1:end-1) ./ h_list(2:end));
fprintf('%.2f  ', space_error_rates);

% plot
figure;
loglog(h_list, spatial_errors, 'bo-', 'LineWidth', 2);
hold on;
loglog(h_list, spatial_errors(1)*(h_list/h_list(1)).^2, 'r--');
xlabel('h');
ylabel('Error');
title('Space Convergence: O(h^2)');
legend('Numerical', 'O(h^2)', 'Location', 'best');
grid on;

%% initial fixed spatial step settings
n_fixed = 500;
h_fixed = 1 / (n_fixed+1);
x_fixed = (1:n_fixed)' * h_fixed;

% time step setting
T = 0.05;
tss_list = [0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0014625];

% initial condition
phi_fixed = sin(2*pi*x_fixed);


% pre-set
time_errors = zeros(size(tss_list));

% build matrix
D_fixed = (diag(ones(n_fixed-1,1),1) + diag(-2*ones(n_fixed,1)) + diag(ones(n_fixed-1,1),-1)) * (1/h_fixed^2);

for i = 1:length(tss_list)
    tss = tss_list(i);
    nsteps = ceil(T/tss); % Ensure that the number of steps is an integer

    % Using Crank Nicholson
    I = eye(n_fixed);
    B = I - (tss/2)*D_fixed;
    L = I + (tss/2)*D_fixed;
    U = phi_fixed;

    % calcu the numerical solution for internal points at time t
    for step = 1:nsteps
        U = B \ (L*U);
    end

    % excat soluction
    U_exact = exp(-4*pi^2*T) .* sin(2*pi*x_fixed);

    % time error calcu and record
    time_errors(i) = max(abs(U - U_exact));
end

% calcu the time convergence rates
time_error_rates = log(time_errors(1:end-1) ./ time_errors(2:end)) ./ log(tss_list(1:end-1) ./ tss_list(2:end));
fprintf('%.2f  ', time_error_rates);

% plot
figure;
loglog(tss_list, time_errors, 'bo-', 'LineWidth', 2);
hold on;
loglog(tss_list, time_errors(1)*(tss_list/tss_list(1)).^2, 'r--');
xlabel('tss');
ylabel('Error');
title('Time Convergence: O(h^2)');
legend('Numerical', 'O(h^2)', 'Location', 'best');
grid on;

%% 
% To verify convergence on spatial and time, I separately fixed the time
% step size and spatial step size. The final results show that the algorithm converges
% both on space and time in $O(h^2)$ time.