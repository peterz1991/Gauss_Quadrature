function [J, D] = symmetrize_T(T)
%-------------------------------------
% Symmetrize a tri-diagnoal matrix
% J = D^{-1}TD, D is a diagonal matrix
%
%
% Open source code, provided by 
% Yue Zhang, 12/03/2016.
% Department of Mathematics, Applied Mathematics and Statistics.
% Case Western Reserve University
%-------------------------------------
a = diag(T);
w = diag(T,1);
b = diag(T,-1);
d = zeros(size(a));

% Calculate the diagonal entries of D
d(2:end) = cumprod(b)./cumprod(w);
d(2:end) = sqrt(d(2:end));
d(1) = 1;
D = diag(d);
Dinv = diag(1./d);
% Calculate the symmetric Jacobi matrix
J = Dinv*T*D;
