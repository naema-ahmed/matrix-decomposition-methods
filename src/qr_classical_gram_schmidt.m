function [Q, R] = hw2_qr(A)
% STUDENT ID: 100064664
% QR Factorisation using Classical Gram-Schmidt (CGS)
% INPUT:
%     A: (m x n) double
% OUTPUT:
%     Q: (m x n) with orthogonal columns
%     R: (n x n) upper-triangular matrix

% Initialisations
[m, n] = size(A);
Q = zeros(m, n);
R = zeros(n, n);
R(1, 1) = sqrt(sum(A(:, 1).*A(:, 1))) ; % R(1,1) is the norm of A in CGS
Q(:, 1) =  A(:, 1)/R(1, 1) ; % first vector of Q is normalised first vector of A

% Looping through the rest of the columns to compute Q (CGS)
for i = 2:n 
    summ = zeros(m,1); % the sum of the projections
    for j = 1:i-1
        summ = summ + (A(:, i)'*Q(:, j))*Q(:, j); % the projection of A column onto the previous Q column
    end
    Q(:, i) = A(:, i) - summ; % subtracting the projections to form the new (i-th) orthogonal column
    Q(:, i) = Q(:, i)/sqrt(sum(Q(:, i).*Q(:, i))); % normalise to make it orthonormal
end

% Computing R
for i = 1:n
    for j = i:n
        R(i, j) = A(:, j)'*Q(:, i); % projections of A columns onto the Q columns
    end 
end 

end
