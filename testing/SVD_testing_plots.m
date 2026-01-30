% SVD Plots
% STUDENT ID: 100064664

n = [8, 64, 256]; % 512 taking too long; code fails to end

for i = 1:1:length(n)
    A = randn(n(i));
    tic; % starts timer
    [U,s,Vt] = hw2_svd(A);
    time(i)=toc; % stops timer & records it
    S = diag(s);
    error(i) = norm(A-(U*S*Vt),'fro')/norm(A, 'fro'); % calculating relative constructive error (actual - prediction)/actual
    U_orthog_error(i) = norm((U')*U - eye(n(i)), 'fro');
    V_orthog_error(i) = norm((Vt')*Vt - eye(n(i)), 'fro');
end 

%% time vs size plot
plot(n, time);
xlabel('Size of Matrix (n)');
ylabel('Time Taken (seconds)');
title('SVD Time vs. Size');

%% error vs size plot
semilogy(n, error);
xlabel('Size of Matrix (n)');
ylabel('Relative Constructive Error');
title('SVD Relative Constructive Error vs. Size');

%% U orthogonality error vs size plot
semilogy(n, U_orthog_error);
xlabel('Size of Matrix (n)');
ylabel('Orthogonality Error ||U^TU - I||_F');
title('SVD U Orthogonality Error vs. Size');

%% V orthogonality error vs size plot
semilogy(n, V_orthog_error);
xlabel('Size of Matrix (n)');
ylabel('Orthogonality Error ||V^TV - I||_F');
title('SVD V Orthogonality Error vs. Size');

%% singular val decay plot

A_tall = randn(100,10);
A_fat = randn(10,100);
[~, s_tall, ~] = hw2_svd(A_tall);
[~, s_fat, ~] = hw2_svd(A_fat);

semilogy(s_tall, 'o-'); hold on;
semilogy(s_fat, 'x-');
xlabel('Index (k)');
ylabel('Singular value (σ_k)');
legend('Tall (100×10)','Wide (10×100)');
title('Singular Value Decay for Short-Fat (10x100) Matrix');
hold off;
 

%% Best Rank-k Approximation Error
shapes = {[100,10], [10,100]};    % tall-skinny and short-fat shapes

shape = 1; % tall-skinny case
    m = shapes{shape}(1); n = shapes{shape}(2);
    A = randn(m,n);
    [U,s,Vt] = hw2_svd(A);
    S = diag(s);
    errors = zeros(length(s),1);
    for k = 1:length(s)
        Ak = U(:,1:k)*S(1:k,1:k)*Vt(1:k,:);
        errors(k) = norm(A - Ak,'fro');
    end
   
    semilogy(1:length(s), errors);
    xlabel('Rank k'); ylabel('||A - A_k||_F');
    title(sprintf('Best Rank-k Approximation Error (%dx%d)',m,n));
    
shape=2; % short-fat case
 m = shapes{shape}(1); n = shapes{shape}(2);
    A = randn(m,n);
    [U,s,Vt] = hw2_svd(A);
    S = diag(s);
    errors = zeros(length(s),1);
    for k = 1:length(s)
        Ak = U(:,1:k)*S(1:k,1:k)*Vt(1:k,:);
        errors(k) = norm(A - Ak,'fro');
    end
    
    semilogy(1:length(s), errors);
    xlabel('Rank k'); ylabel('||A - A_k||_F');
    title(sprintf('Best Rank-k Approximation Error (%dx%d)',m,n));