 function [U,s,Vt] = hw2_svd(A)
 % STUDENT ID: 100064664
 % Single Value Decomposition  via eigen-route on A'*A using my Schur
 %  INPUT:
 % A : (m x n) double
 % OUTPUT:
 % U : (m x r) orthonormal columns (economy)
 % s : (r x 1) non-negative singular values (sorted)
 % Vt: (r x n) 

 % Initialisations & using hw2_schur
 [m, n] = size(A);
 [V, T] = hw2_schur(A'*A, 5000, 1e-10, false); % Since AtA is symmetric, shifting is not necessary, and leads to near 0 values as pivots
 s = diag(T);
 s(s<0) = 0;
 s = sqrt(s);
 

 % Sort s in decsending order & reorder Vt to match corresponding singular values
 [s, index] = sort(s,"descend");
 Vr = V(:, index);
 Vt = Vr';
 
 % Trunkating s for singular values smaller than an arbitrary threshold
 k = max(size(A)) * eps(norm(A)); % Define a threshold for small singular values
 s = s(s >= k);
 r = length(s);

 % Trunkate Vt to be (r x n)
 Vt = Vt(1:r,:);


 % Instead of finding full U from U=AVs, find U for only large enough s
 U = zeros(m, length(s));
 for i = 1:length(s)
     if s(i) > k
         U(:,i) = A * V(:, i) * (1/s(i));
     end
 end

 end
