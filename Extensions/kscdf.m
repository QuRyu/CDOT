% © Brent R. Hickman and Timothy P. Hubbard
% hickmanbr@uchicago.edu and timothy.hubbard@colby.edu

% please cite our Journal of Applied Econometrics paper "Replacing Sample
% Trimming with Boundary Correction in Nonparametric Estimation of
% First-Price Auctions" if this code is helpful to you.

% construct kernel cdf or edf for an observed vector

function [F,evalgrid] = kscdf(obs,kernel,evalgrid)

n = length(obs);

if (nargin == 1)
	kernel = 'edf';
end

% points to evaluate edf at (if not provided)
if (nargin == 2)
    evalgrid = linspace(min(obs) - 1,max(obs) + 1,n*10);
end

% use Hardle's (1991) bandwidth transformation constant
c = 1;
if (isequal(kernel,'gaussian') == 1)
    c = 1;
elseif (isequal(kernel,'epanechnikov') == 1)
    c = 2.214;
elseif strcmpi(kernel,'triweight')
    c = 2.978;
elseif (isequal(kernel,'uniform') == 1)
    c = 1.740;
elseif (isequal(kernel,'triangle') == 1)
    c = 2.432;
elseif (isequal(kernel,'cosinus') == 1)
    c = 2.288;
elseif (isequal(kernel,'quartic') == 1)
    c = 2.623;
end

% use Silverman's (1986) optimal bandwidth
h = c*(4/3)^(1/5)*std(obs)*(n)^(-1/5);

% create kernel cdfs or edfs
T = length(evalgrid);
F = zeros(T, 1);
for t=1:T
	u = (evalgrid(t) - obs)/h;
	if strcmpi(kernel,'edf')
		F(t) = length(find(evalgrid(t) >= obs));
	end
end
F = F/n;
