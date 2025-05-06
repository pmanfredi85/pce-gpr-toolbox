function P = orthonormal_hermite(k,N)

% P = orthonormal_hermite(k,N) returns the coefficients of the k-th
% normalized Hermite polynomial. N is optional and indicates the maximum
% degree of the output polynomial. Of course, N must be >= k. If not
% specified, N = k is taken. Otherwise, the polynomial P is zero padded
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

% Check input arguments
if nargin==1
    N = k;
elseif nargin==2
    if N<k
        error('N must be equal or greater than k');
    end
else
    error('Wrong number of input arguments');
end

% H_(k-2)
Pm2 = zeros(N+1,1);
Pm2(end) = 1;

if k==0
    P = Pm2;
    return;
end

% H_(k-1)
Pm1 = zeros(N+1,1);
Pm1(end-1) = 1;
    
if k==1
    P = Pm1;
    return;
end

% Iterative calculation of H_k though recursion relation
y = zeros(N+1,1);
y(end-1) = 1;
    
for I = 2:k
    Paux = conv(y,Pm1);
    P = Paux(end-N:end) - (I-1)*Pm2;
    Pm2 = Pm1;
    Pm1 = P;
end

% Orthonormalization
P = P/sqrt(factorial(k));