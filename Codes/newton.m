function [x, res, xvec, resvec] = newton(f, df, x0, maxiter, tol)
% NEWTON  Newton's method for solving f(x) = 0
%
%   Usage:
%       [x, res, xvec, resvec] = newton(f, df, x0, maxiter, tol)
%
%   Inputs:
%       f        - function handle for f(x)
%       df       - function handle for f'(x)
%       x0       - initial guess (scalar)
%       maxiter  - maximum number of iterations allowed (positive integer)
%       tol      - tolerance for successive difference |x_{k+1}-x_k| (positive scalar)
%
%   Outputs:
%       x        - final approximation to a root x*
%       res      - residual |f(x)| at the final approximation
%       xvec     - column vector storing the iterates {x_k}, including x0 and the final x
%       resvec   - column vector storing {|f(x_k)|}, same length as xvec

    % ---------- Input checking ----------
    % Validate the number of inputs: the function requires exactly 5 inputs.
    if nargin ~= 5
        error('newton:InvalidNargin', 'Exactly five inputs are required: f, df, x0, maxiter, tol.');
    end
    
     % Check that both f and df are valid function handles.
    if ~isa(f, 'function_handle') || ~isa(df, 'function_handle')
        error('newton:InvalidHandle', 'f and df must be valid function handles.');
    end

    % Ensure the initial guess x0 is a finite scalar value.
    if ~isscalar(x0) || ~isfinite(x0)
        error('newton:InvalidX0', 'x0 must be a finite scalar.');
    end

    % Verify that maxiter is a positive integer.
    if ~isscalar(maxiter) || ~isfinite(maxiter) || maxiter <= 0 || floor(maxiter) ~= maxiter
        error('newton:InvalidMaxIter', 'maxiter must be a positive integer.');
    end

    % Verify that the tolerance tol is a positive, finite scalar.
    if ~isscalar(tol) || ~isfinite(tol) || tol <= 0
        error('newton:InvalidTol', 'tol must be a positive finite scalar.');
    end

    % ---------- Preallocation & init ----------
    xvec   = zeros(maxiter+1, 1);   % store x0 ... up to at most x_maxiter
    resvec = zeros(maxiter+1, 1);   % store |f(x_k)|
    k = 0;
    xk     = x0;
    xvec(1)   = xk;
    resvec(1) = abs(f(xk));
    converged = false;

    % Print header
    fprintf(' Newton''s method iterations\n');
    fprintf('    k\t\t   x_k\t          |f(x_k)|      |x_{k+1}-x_k|\n');
    fprintf(' ---------------------------------------------------------------------------\n');

    % ---------- Main iteration loop ----------
    while k < maxiter
        fk  = f(xk);
        dfk = df(xk);

        % Check derivative is not (near-)zero
        if ~isfinite(dfk) || abs(dfk) < eps
            fprintf('Derivative is zero or ill-conditioned at k=%d (x=%.16g). Stopping.\n', k, xk);
            break;
        end

        % Newton step update
        xnew = xk - fk/dfk;
        step = abs(xnew - xk);

        % Print current row 
        fprintf(' %4d\t% .16e\t% .3e\t% .3e\n', k, xk, abs(fk), step);

        % Accept new iterate
        k = k + 1;
        xk = xnew;

        xvec(k+1)   = xk;  % store new iterate
        resvec(k+1) = abs(f(xk));

        % Check stopping criterion on successive difference
        if step < tol
            converged = true;
            break;
        end
    end

    % ---------- Trim output vectors to actual length ----------
    xvec   = xvec(1:k+1);
    resvec = resvec(1:k+1);

    % ---------- Final outputs ----------
    x   = xk;
    res = abs(f(x));

    % ---------- Final message ----------
    if converged
        fprintf('Convergence achieved: |x_{k+1}-x_k| < tol after %d iterations. Final |f(x)| = %.3e\n', k, res);
    else
        if k >= maxiter
            fprintf('Maximum iterations reached (%d). Final |f(x)| = %.3e\n', maxiter, res);
        else
            fprintf('Terminated early due to derivative breakdown. Final |f(x)| = %.3e\n', res);
        end
    end
end