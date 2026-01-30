% Schur Plots
% STUDENT ID: 100064664

n = [8, 64, 256]; % 512 taking too long; code fails to end

for i = 1:1:length(n)
    A = randn(n(i));
    tic; % starts timer
    [Q, T] = hw2_schur(A);
    time(i)=toc; % stops timer & records it
    error(i) = norm(A-(Q*T*(Q')),'fro')/norm(A, 'fro'); % calculating relative constructive error (actual - prediction)/actual
    orthog_error(i) = norm((Q')*Q - eye(n(i)), 'fro');
end
%%
clear; clc;
n=64;
A=randn(n);
tol = 1e-10;
max_iter = 5000;
use_shift = true;
Qn = eye(n);
belowdiag = tril(A,-1);
toler = norm(belowdiag, 'fro'); % Using norm of elements below diagonal to measure convergence
iter=0;
tolers=[];
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
    tolers = [tolers toler];
    iter = iter + 1;
end 

%% time vs size plot
plot(n, time);
xlabel('Size of Matrix (n)');
ylabel('Time Taken (seconds)');
title('Schur Time vs. Size');

%% error vs size plot
semilogy(n, error);
xlabel('Size of Matrix (n)');
ylabel('Relative Constructive Error');
title('Schur Relative Constructive Error vs. Size');

%% orthogonality error vs size plot
semilogy(n, orthog_error);
xlabel('Size of Matrix (n)');
ylabel('Orthogonality Error ||Q^TQ - I||_F');
title('Schur Orthogonality Error vs. Size');

%% convergence: subdiag norm vs iterations
semilogy(tolers);
xlabel('Iteration');
ylabel('||subdiagonal A_k||_F');
title('Schur by QR Convergence (n=64)')