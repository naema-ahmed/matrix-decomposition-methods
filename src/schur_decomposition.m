 function [Q,T] = hw2_schur(A, max_iter, tol, use_shift)
 % STUDENT ID: 100064664
 % Real Schur Factorisation via QR Algorithm
 % INPUT:
 %     A: (n x n) double (real)
 %     max_iter: integer (default 5000)
 %     tol: double (default 1e-10)
 %     use_shift: logical (default true)
 % OUTPUT:
 %     Q : (n x n) orthogonal
 %     T : (n x n) real (quasi-)upper-triangular (Schur form)
 % NOTES:
 %     QR iterations; symmetric A is acceptable as a starting scope

 %  Assigning default values to missing arguments

 if nargin == 3
     use_shift = true;
 elseif nargin == 2
     tol = 1e-10;
     use_shift = true;
 elseif nargin == 1
     max_iter = 5000;
     tol = 1e-10;
     use_shift = true;
 end 

% Initialisations

[~,n] = size(A); 
Qn = eye(n);
belowdiag = tril(A,-1);
toler = norm(belowdiag, 'fro'); % Using norm of elements below diagonal to measure convergence
iter=0;

% Loop to perform Schur Decomposition


while (toler >= tol) && (iter <= max_iter)

    if use_shift % When true, use shift by bottom right element (Rayleigh)
        s = A(n,n);
        if isfinite(s)==false % To prevent NaN or +/- Inf caused by my QR, which may divide by 0 or near 0
            s = 0;
        end 

    else % When false, no shifting (same as shift by 0)
        s=0; 
    end

    A = A - s*eye(n); % Shift A matrix
    [Qk, Rk] = hw2_qr(A); % Perform my CGS QR on shifted A
    A = Rk*Qk + s*eye(n); % Use similarity and undo the shift 
    Qn = Qn*Qk; % Record total product of Qk rotations in Qn
    belowdiag = tril(A,-1); % Extract elements below diagonal of A
    toler = norm(belowdiag, 'fro'); % Update toler (shows 'error' since I want below diagonal to be near 0 for quasi-upper triangular)
    iter = iter + 1;
end 

% In case max_iter reached before toler <= tol

if iter == max_iter
    fprintf('Maximum iterations reached before convergence. ( toler at max_iter = %.3e )\n', toler)
end 

    Q = Qn;
    T = A;

 end 


