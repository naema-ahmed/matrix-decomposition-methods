function [L,U] = hw2_lu(A)
% STUDENT ID: 100064664
% LU Factorisation by Doolittle's Method
% INPUT:
%     A: (n x n)
% OUTPUT:
%     L: (n x n) unit lower-triangular matrix (diag=1)
%     U: (n x n) upper-triangular matrix

[~, n] = size(A);
L = eye(n);
U = zeros(n);

% Initialise first row of U & first column of L
U(1, :) = A(1, :);
L(:, 1) = A(:, 1)/U(1, 1);

% Doolittle algorithm loop for all other elements
for i = 2:n
    for j = i:n
        U(i, j) = A(i, j) - L(i, 1:i-1) * U(1:i-1, j);
    end 

    if abs(U(i,i)) < eps % Detecting zero/near-zero pivots
        error("Zero or near-zero pivot detected at iteration %d", i)
    end 

    for k=i+1:n
        L(k, i) = ( A(k, i) - L(k, 1:i-1) * U(1:i-1, i) )/ U(i, i);
    end
end 
end 
