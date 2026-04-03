function [tHist, uHist] = ForwardEuler(f, tspan, u0, nsteps)
% FORWARDEULER  Forward Euler method for solving ODEs
% y'(t) = f(t, y(t)) for t in [t0, T] with initial condition y(t0)=y0.
% Inputs:
%       f: a function handle f for the function f of (3);
%       tspan: a 2-vector tspan, with components t0 and T;
%       u0: a p × 1 vector u0, containing the initial condition y0;
%       nsteps: the number of steps nteps to be taken;
% Outputs:
%       tHist: an (m + 1) × 1 vector tHist, storing the approximation times ti
%       uHist: an (m + 1) × p matrix uHist, whose nth row stores the approximation ui ≈ y(ti) ∈ Rp
    
    % check input
    if nargin < 4
        error('Not enough input arguments');
    end
    if length(tspan) ~= 2
        error('tspan must be a 2-vector');
    end
    if nsteps <= 0
        error('nsteps must be positive');
    end
   

    % extract time info
    t0 = tspan(1);
    T = tspan(2);

    % time step size
    tss = (T-t0)/nsteps;

    % pre-set memory
    p = length(u0);
    tHist = zeros(nsteps+1,1);
    uHist = zeros(nsteps+1,p);

    % initial
    tHist(1) = t0;
    uHist(1,:) = u0';

    % time step loop
    for i = 1:nsteps
        t_current = tHist(i);
        u_current = uHist(i,:)';

        % update
        u_new = u_current + tss * f(t_current, u_current);
        tHist(i+1) = t_current + tss;
        uHist(i+1,:) = u_new';
    end
end
