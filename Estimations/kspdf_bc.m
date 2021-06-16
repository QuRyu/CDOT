% � Brent R. Hickman and Timothy P. Hubbard
% hickmanbr@uchicago.edu and timothy.hubbard@colby.edu

% please cite our Journal of Applied Econometrics paper "Replacing Sample
% Trimming with Boundary Correction in Nonparametric Estimation of
% First-Price Auctions" if this code is helpful to you.

% compute boundary-corrected kernel smoothed estimate of pdf

function [fhat,h] = kspdf_bc(obs,kernel,evalgrid,b0,usehardle)

% number of points in observed data
n = length(obs);

% use Hardle's (1991) bandwidth transformation constant
c = 1;
if (isequal(kernel,'gaussian') == 1)
    c = 1;
elseif (isequal(kernel,'epanechnikov') == 1)
    c = 2.214;
elseif (isequal(kernel,'triweight') == 1)
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
if (exist('usehardle','var'))
	if usehardle == 0
		c = 1;
	elseif usehardle == 1
		c = c;
	end
end
h = c*(4/3)^(1/5)*std(obs)*(n)^(-1/5);

% shift support
lowobs = min(obs);
obs = obs - lowobs;
evalgrid = evalgrid - lowobs;

% h1 bandwidth parameter used in dhat and everything underneath it---see
% (2.9) in KZ, line that follows it, and (2.3), (2.4), and (2.5) which are
% explicit below 
h1 = h*n^(-1/20);

% constant A must be such that 3A > 1.  A = .55 was suggested for practical
% use by KZ pg 503 and equation (2.9)
A = 0.55;
if nargin == 3
	% b0 is needed for h0 bandwidth---see KZ (2008) bottom of page 501 and
	% the K0 term in (4.9) of ZK (1998) which is K0 = 6 + 18t + 12t^2,
	% t\in[-1,0]; this will compute b0 for any kernel choice provided by
	% user 
	num_t1 = quad(@evalkspdf_num,-1,1,1e-14,[],kernel);
	num_t2 = 0 - (-144/5 + 432/4 - 468/3 + 216/2 - 36);
	denom_t1 = 0 - (-2 + 18/4 - 12/5);
	denom_t2 = quad(@evalkspdf_denom,-1,1,1e-14,[],kernel);
	b0 = ((num_t1^2*num_t2)/(denom_t1^2*denom_t2))^(1/5);
end
% third bandwidth parameter for estimation of dhat, see KZ pg 501
h0 = b0*h1;

% evaluate fstar_h1---see (2.4) of KZ but use h1 instad of h
u1 = (h1 - obs)/h1;
k1 = evalkspdf(u1,kernel);
fstar_h1 = sum(k1)/(n*h1);
f_h1 = fstar_h1 + 1/(n^2);

% evaluate fstar_0---see (2.5) of KZ and (4.9) of ZK
u0 = -obs/h0;
k0 = (6 + 18*u0 + 12*u0.^2).*(u0 >= -1).*(u0 <= 0);
fstar_0 = sum(k0)/(n*h0);
f_0 = max(fstar_0,1/(n^2));

% evaluate dhat---see (2.3) of KZ and line below (2.9) of KZ
dhat = (log(f_h1) - log(f_0))/h1;

% evaluate observed data points (or pseudo values) at ghat---see (2.9) of
% KZ
ghat_obs = obs + dhat*obs.^2 + A*dhat^2*obs.^3;

T = length(evalgrid);
% for density = 0 rule beyond min and max observations
maxobs = max(obs);
minobs = min(obs);
fhat = zeros(T, 1);

for t=1:T
	u = (evalgrid(t) - obs)/h;
	ku = evalkspdf(u,kernel);
	ug = (evalgrid(t) + ghat_obs)/h;
    kug = evalkspdf(ug,kernel);
    kterm = ku + kug;
	fhat(t) = sum(kterm);
    % for density = 0 rule beyond min and max observations
    if (evalgrid(t) > maxobs) || (evalgrid(t) < minobs)
        fhat(t) = 0;
    end
end
fhat = fhat/(n*h);

% shift support back---doesn't matter now but if we ever output these
% values it would
obs = obs + lowobs;
evalgrid = evalgrid + lowobs;


