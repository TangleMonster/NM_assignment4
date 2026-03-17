function [x, res, xvec, resvec] = secant(f, x0, x1, maxiter, tol)
% SECANT  Secant method for solving f(x) = 0
%
%   Usage:
%       [x, res, xvec, resvec] = secant(f, x0, x1, maxiter, tol)
%
%   Inputs:
%       f        - function handle for f(x)
%       x0       - first initial guess (scalar)
%       x1       - second initial guess (scalar)
%       maxiter  - maximum number of iterations allowed (positive integer)
%       tol      - tolerance for successive difference |x_{k+1}-x_k| (positive scalar)
%
%   Outputs:
%       x        - final approximation to a root x*
%       res      - residual |f(x)| at the final approximation
%       xvec     - column vector storing the iterates {x_k}, including x0 and the final x
%       resvec   - column vector storing {|f(x_k)|}, same length as xvec

    % ---------- Input checking ----------
    % Validate the number of inputs: exactly 5 are required.
    if nargin ~= 5
        error('secant:InvalidNargin', ...
              'Exactly five inputs are required: f, x0, x1, maxiter, tol.');
    end

    % Check that f is a valid function handle.
    if ~isa(f, 'function_handle')
        error('secant:InvalidHandle', 'f must be a valid function handle.');
    end

    % Ensure initial guesses are finite scalars.
    if ~isscalar(x0) || ~isfinite(x0)
        error('secant:InvalidX0', 'x0 must be a finite scalar.');
    end
    if ~isscalar(x1) || ~isfinite(x1)
        error('secant:InvalidX1', 'x1 must be a finite scalar.');
    end

    % Verify maxiter is a positive integer.
    if ~isscalar(maxiter) || ~isfinite(maxiter) || maxiter <= 0 || floor(maxiter) ~= maxiter
        error('secant:InvalidMaxIter', 'maxiter must be a positive integer.');
    end

    % Verify tol is a positive, finite scalar.
    if ~isscalar(tol) || ~isfinite(tol) || tol <= 0
        error('secant:InvalidTol', 'tol must be a positive finite scalar.');
    end

    % ---------- Preallocation & init ----------
    % Allocate for up to maxiter+2 entries
    xvec   = zeros(maxiter+2, 1);
    resvec = zeros(maxiter+2, 1);
    xkm1 = x0;     fkm1 = f(xkm1);   % k-1 state
    xk   = x1;     fk   = f(xk);     % k state
    xvec(1)   = xkm1;  resvec(1) = abs(fkm1);
    xvec(2)   = xk;    resvec(2) = abs(fk);
    converged = false;
    k = 1;  % corresponds to x1 being the "current" iterate

    % Print header (k counts the "current" iterate index)
    fprintf(' Secant method iterations\n');
    fprintf('    k\t\t   x_k\t          |f(x_k)|      |x_{k+1}-x_k|\n');
    fprintf(' ---------------------------------------------------------------------------\n');

    % ---------- Main iteration loop ----------
    while k <= maxiter
        denom = (fk - fkm1);

        % Check derivative is not (near-)zero
        if ~isfinite(denom) || abs(denom) < eps
            fprintf('Secant slope is zero/ill-conditioned at k=%d (x_k=%.16g). Stopping.\n', k, xk);
            break;
        end

        % Secant update
        xnew = xk - fk * (xk - xkm1) / denom;
        step = abs(xnew - xk);

        % Print current row
        fprintf(' %4d\t% .16e\t% .3e\t% .3e\n', k, xk, abs(fk), step);

        % Accept new iterate
        k      = k + 1;
        xkm1   = xk;     fkm1 = fk;
        xk     = xnew;   fk   = f(xk);

        xvec(k+1)   = xk;  % store new iterate
        resvec(k+1) = abs(fk);

        % Check stopping criterion
        if step < tol
            converged = true;
            break;
        end
    end

    % ---------- Trim output vectors ----------
    % We stored up to index (k+1); remove trailing zeros.
    last_idx = find(xvec ~= 0 | resvec ~= 0, 1, 'last');
    if isempty(last_idx), last_idx = 2; end
    xvec   = xvec(1:last_idx);
    resvec = resvec(1:last_idx);

    % ---------- Final outputs ----------
    x   = xk;
    res = abs(f(x));

    % ---------- Final message ----------
    if converged
        fprintf('Convergence achieved: |x_{k+1}-x_k| < tol after %d iterations. Final |f(x)| = %.3e\n', k-1, res);
    else
        if k > maxiter
            fprintf('Maximum iterations reached (%d). Final |f(x)| = %.3e\n', maxiter, res);
        else
            fprintf('Terminated early due to ill-conditioned secant slope. Final |f(x)| = %.3e\n', res);
        end
    end
end
