function Jk = Jacobi(k,polynomial_name)
%----------------------------------------
% Construct a k x k Jacobi matrix of certain ONPs.
%----------------------------------------

if strcmp(polynomial_name,'Legendre')==1
    % construct the Jacobi matrix of Legendre polynomials of order n
    d = 0:k-1;
    d1 = (d+1)./(2*d+1);    % Construct the upper diag entries of Jk
    d2 = d./(2*d+1);        % Construct the lower diag entries of Jk
    Tk = diag(d1(1:end-1),1) + diag(d2(2:end),-1);
    Jk = symmetrize_T(Tk);
    
elseif strcmp(polynomial_name,'Laguerre')==1
    d = 0:k-1;
    d1 = -(d+1);
    d2 = 2*d+1;
    d3 = -d;
    Jk = diag(d1(1:end-1),1) + diag(d3(2:end),-1) + diag(d2);
    
elseif strcmp(polynomial_name,'Chebyshev1')==1
    d = 1/2*ones(k-1,1);
    Jk = diag(d,1) + diag(d,-1);
    if k>=2
        Jk(1,2) = 1;
    end
    Jk = symmetrize_T(Jk);

elseif strcmp(polynomial_name,'Chebyshev2')==1
    d = 1/2*ones(k-1,1);
    Jk = diag(d,1) + diag(d,-1);
end