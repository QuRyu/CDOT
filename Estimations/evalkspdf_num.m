% © Brent R. Hickman and Timothy P. Hubbard
% hickmanbr@uchicago.edu and timothy.hubbard@colby.edu

% please cite our Journal of Applied Econometrics paper "Replacing Sample
% Trimming with Boundary Correction in Nonparametric Estimation of
% First-Price Auctions" if this code is helpful to you.

% evaluate specified kernel density at a grid of points

function k = evalkspdf(u,kernel)

if (isequal(kernel,'gaussian') == 1)
	k = 1/sqrt(2*pi)*exp(-(u.^2)/2);
elseif (isequal(kernel,'epanechnikov') == 1)
	k = (3/4*(1 - u.^2)).*(abs(u) <= 1);
elseif (isequal(kernel,'triweight') == 1)
	k = (35/32*(1 - u.^2).^3).*(abs(u) <= 1);
elseif (isequal(kernel,'uniform') == 1)
	k = (ones(size(u))*1/2).*(abs(u) <= 1);
elseif (isequal(kernel,'triangle') == 1)
	k = (1 - abs(u)).*(abs(u) <= 1);
elseif (isequal(kernel,'cosinus') == 1)
	k = (pi/4*cos(pi/2*u)).*(abs(u) <= 1);
elseif (isequal(kernel,'quartic') == 1)
	k = (15/16*(1 - u.^2).^2).*(abs(u) <= 1);
else
	fprintf('Kernel choice not defined or spelt wrong\n')
	k = -Inf;
end

% for computing b0
k = u.^2.*k;