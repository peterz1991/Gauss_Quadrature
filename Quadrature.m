function I = Quadrature(f,n,type) 
%------------------------------------------
% A comprehensive quadrature rule list with n points for function f.
% We assume the integration interval is [-1,1].
% type is a string variable controling which rule to implement:
% type = 'Gauss'   // Gauss rules with Legendre polynomials;
%      = 'Laguerre' // rules with Laguerre polynomials;
%      = 'Chebyshev1' // rules with Chebyshev polynomials of 1st kind;
%      = 'Chebyshev2' // rules with Chebyshev polynomials of 2nd kind;
%      = 'Gauss_Radau' // Gauss_Radau rules (one more node on boundary);
%      = 'Gauss_Lobatto' // Gauss_Lobatto rules (two more nodes);
%      = 'Anti_Gauss' // Anti-Gauss Rules;
%
%
% Open source code, provided by 
% Yue Zhang, 12/03/2016.
% Department of Mathematics, Applied Mathematics and Statistics.
% Case Western Reserve University
%-------------------------------------------

if strcmp(type,'Gauss') == 1
    % n-pt Gauss quadrature of f (with Legendre polynomials)
    % modified from [1] (see references)
    beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
    J = diag(beta,1) + diag(beta,-1); % Jacobi matrix
    % An alternative way to obtain J
    % J = Jacobi(n,'Legendre');
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes (= Legendre points)
    w = 2*V(1,i).^2; % weights 
    I = w*f(x); % the integral from -1 to 1.

elseif strcmp(type,'Laguerre') == 1
    J = Jacobi(n,'Laguerre');
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes (= Legendre points)
    w = V(1,i).^2; % Laguerre weight
    g = @(x) f(x).*exp(x); % modify f with corrspd weights
    I = w*g(x);  
    
elseif strcmp(type,'Chebyshev1') == 1
    J = Jacobi(n,'Chebyshev1');
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes xj = cos((2j-1)*pi/(2n)), j=1...n
    w = pi*V(1,i).^2;  % w = pi/n (constant for all nodes)
    g = @(x) f(x).*sqrt(1-x.^2); % modify f with corrspd weights
    I = w*g(x);  
    % Note*: for Chebyshev1&2, explicit sol for x and w exists.

elseif strcmp(type,'Chebyshev2') == 1
    J = Jacobi(n,'Chebyshev2');
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes xj = cos(j*pi/(n+1)), j=1...n 
    w = pi/2*V(1,i).^2; % wj = pi/(n+1)*sin^2(j*pi/(n+1)), j=1...n
    g = @(x) f(x)./sqrt(1-x.^2); % modify f with corrspd weights
    I = w*g(x);  
    
elseif strcmp(type,'Gauss_Radau') == 1
    % n-pt Gauss-Radau rule + 1 fixed node at -1 (with Legendre polys)    
    J = Jacobi(n+1,'Legendre');
    a = 1; % fixed node
    beta_N = J(n,n+1);
    A = J(1:n,1:n) - a*eye(n);
    b = [zeros(n-1,1);beta_N^2];
    z = A\b;
    alpha_nplus1 = z(end) + a; % construct alpha_{n+1} = a + delta(end)
    J(n+1,n+1) = alpha_nplus1; % Extended Jacobi matrix
    
    % proceed to gauss rule as before
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes (= Legendre points)
    w = 2*V(1,i).^2; % generalized Laguerre weight
    I = w*f(x);   
    
elseif strcmp(type, 'Gauss_Lobatto') == 1
    % n-pt Gauss-Lobatto rule + two fixed nodes at +-1 (with Legendre polys)
    a = -1; b = 1; % fixed nodes
    beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
    J = diag(beta,1) + diag(beta,-1); % Jacobi matrix
    I = eye(n);
    eN = I(:,end);
    gamma_vector = (J-a*eye(n))\eN;
    mu_vector = (J-b*eye(n))\eN;
    A = [1 -gamma_vector(end); 1 -mu_vector(end)];
    b = [a; b];
    z = A\b;    % returns a 2 x 1 vector with extended weights
    J_extend = zeros(n+1);  % start to construct extended Jacobi
    J_extend(1:n,1:n) = J;
    J_extend(n+1,n+1) = z(1);
    J_extend(n,n+1) = sqrt(z(2));
    J_extend(n+1,n) = sqrt(z(2));
    J = J_extend;
    
    % proceed to gauss rule as before
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes (= Legendre points)
    w = 2*V(1,i).^2; 
    I = w*f(x);   
    
elseif strcmp(type,'Anti_Gauss') == 1
    % (n+1)-pt Anti-Gauss quadrature rule (with Legendre polynomials)
    J = Jacobi(n+1,'Legendre');
    J(n,n+1) = sqrt(2)*J(n,n+1);
    J(n+1,n) = sqrt(2)*J(n+1,n);
    
    % proceed to gauss rule as before 
    [V,D] = eig(J); % eigenvalue decomposition
    x = diag(D); [x,i] = sort(x); % nodes (= Legendre points)
    w = 2*V(1,i).^2; % weights 
    I = w*f(x); % the integral from -1 to 1.    
    
end   

end
%% References
% [1] Lloyd N. Trefethen, Is Gauss Quadrature Better than Clenshaw–Curtis?.
% [2] Golub & Welsch, Calculation of Gauss Quadrature Rules.