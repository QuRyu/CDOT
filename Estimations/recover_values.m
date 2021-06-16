

function [v, vmin_bc, vmax_bc] = recover_values(b, N) 
bmax = max(b);
bmin = min(b);

usehardle = 1; 

vlow = 0;
vhigh = 1;
neval = 500;
evalpts = linspace(vlow,vhigh,neval);

kernel = 'triweight';

% b0 is needed for h0 bandwidth---see KZ (2008) bottom of page 501 and the
% K0 term in (4.9) of ZK (1998) which is K0 = 6 + 18t + 12t^2, t\in[-1,0];
% this will compute b0 for any kernel choice provided by user
num_t1 = quad(@evalkspdf_num,-1,1,1e-14,[],kernel);
num_t2 = 0 - (-144/5 + 432/4 - 468/3 + 216/2 - 36);
denom_t1 = 0 - (-2 + 18/4 - 12/5);
denom_t2 = quad(@evalkspdf_denom,-1,1,1e-14,[],kernel);
b0 = ((num_t1^2*num_t2)/(denom_t1^2*denom_t2))^(1/5);


[gB_bc,hb_l] = kspdf_bc(b,kernel,b,b0,usehardle);
% right-boundary correction: reflect bids over zero and do left-boundary
% correction 
[gB_bc_right,hb_r] = kspdf_bc(-b,kernel,-b,b0,usehardle);
bind = find(b >= bmax - hb_r);
gB_bc(bind) = gB_bc_right(bind);
GB = kscdf(b,'edf',b);
v_bc = b - (1-GB)./(gB_bc*(N - 1));

v = v_bc';
vmin_bc = min(v_bc);
vmax_bc = max(v_bc);
