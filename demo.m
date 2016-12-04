%
% Open source code, provided by 
% Yue Zhang, 12/03/2016.
% Department of Mathematics, Applied Mathematics and Statistics.
% Case Western Reserve University
%% Functions of test
f1 = @(x) 1./sqrt(1-x.^2); % pi
f2 = @(x) exp(-1./x.^2);   % 2/e - 2*Sqrt(pi)*erfc(1)
f3 = @(x) 1./(1+20*x.^2);  % atan(2*sqrt(5))/sqrt(5)
f4 = @(x) x.^10;      % 2/11
f5 = @(x) exp(x);     % e- 1/e
f6 = @(x) exp(-x.^2); % sqrt(pi)*erf(1)

%% Compute the absolute errors with gauss rule
Nmax = 20; E = zeros(6,Nmax); clf
e = exp(1);

for n = 1:Nmax
    E(1,n) = abs(Quadrature(f1,n,'Gauss') - pi);
    E(2,n) = abs(Quadrature(f2,n,'Gauss') - 2/e + 2*sqrt(pi)*erfc(1));
    E(3,n) = abs(Quadrature(f3,n,'Gauss') - atan(2*sqrt(5))/sqrt(5));
    E(4,n) = abs(Quadrature(f4,n,'Gauss') - 2/11);
    E(5,n) = abs(Quadrature(f5,n,'Gauss') - e + 1/e);
    E(6,n) = abs(Quadrature(f6,n,'Gauss') - sqrt(pi)*erf(1));
end

% Plot results:
labels = {'1/sqrt(1-x^2)','exp(-1/x^{-2})','1/(1+20x^2)'...
    ,'x^{10}','exp(x)', 'exp(-x^2)'};

for iplot = 1:6
    subplot(3,2,iplot)
    semilogy(E(iplot,:)+1e-100,'.','markersize',12), hold on
    plot(E(iplot,:)+1e-100,'linewidth',.8)
    axis([0 Nmax 1e-18 1e3]), grid on
    set(gca,'xtick',0:10:Nmax,'ytick',(10).^(-15:5:0))
    ylabel error, text(32,.0004,labels(iplot))
end