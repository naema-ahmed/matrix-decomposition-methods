% QR Plots
% STUDENT ID: 100064664

n = [8, 64, 256, 512];

for i = 1:1:length(n)
    A = randn(n(i));
    tic; % starts timer
    [Q, R] = hw2_qr(A);
    time(i)=toc; % stops timer & records it
    error(i) = norm(A-(Q*R),'fro')/norm(A, 'fro'); % calculating relative constructive error (actual - prediction)/actual
    orthog_error(i) = norm((Q')*Q - eye(n(i)), 'fro');
end

%% time vs size plot
plot(n, time);
xlabel('Size of Matrix (n)');
ylabel('Time Taken (seconds)');
title('QR Time vs. Size');

%% error vs size plot
semilogy(n, error);
xlabel('Size of Matrix (n)');
ylabel('Relative Constructive Error');
title('QR Relative Constructive Error vs. Size');

%% orthogonality error vs size plot
semilogy(n, orthog_error);
xlabel('Size of Matrix (n)');
ylabel('Orthogonality Error ||Q^TQ - I||_F');
title('QR Orthogonality Error vs. Size');