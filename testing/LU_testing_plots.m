% LU Plots
% STUDENT ID: 100064664

n = [8, 64, 256, 512];

for i = 1:1:length(n)
    A = randn(n(i));
    tic; % starts timer
    [L, U] = hw2_lu(A);
    time(i)=toc; % stops timer & records it
    error(i) = norm(A-(L*U),'fro')/norm(A, 'fro'); % calculating relative constructive error (actual - prediction)/actual
end
%% time vs size plot
plot(n, time);
xlabel('Size of Matrix (n)');
ylabel('Time Taken (seconds)');
title('LU Time vs. Size');

%% error vs size plot
semilogy(n, error);
xlabel('Size of Matrix (n)');
ylabel('Relative Constructive Error');
title('LU Relative Constructive Error vs. Size');