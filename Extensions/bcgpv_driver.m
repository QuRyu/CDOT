% � Brent R. Hickman and Timothy P. Hubbard
% hickmanbr@uchicago.edu and timothy.hubbard@colby.edu

% please cite our Journal of Applied Econometrics paper "Replacing Sample
% Trimming with Boundary Correction in Nonparametric Estimation of
% First-Price Auctions" if this code is helpful to you.

% a sample driver file using uniformly distributed data

clear all
clc

% parameters for figure properties
lwidth = 2;
fsize = 14;
set(0,'defaulttextinterpreter','latex','Defaulttextfontsize',fsize);

% load your data here as a column vector!
% b = data
% N = from your data
% sample uniform data from [0,1/2] with N = 2
% b = rand(600,1);
% b = b/2;
% N = 2;
% bmax = max(b);
% bmin = min(b);
% b = sort(b);

% use hardle bandwidth transformation constant for non-gaussian kernel
usehardle = 1;

% indicator to turn on plotting
plotind = 0;

% indicator to turn on simulation
simulind = 1; 

% choice of kernel---everything else is automated
kernel = 'triweight';

% b0 is needed for h0 bandwidth---see KZ (2008) bottom of page 501 and the
% K0 term in (4.9) of ZK (1998) which is K0 = 6 + 18t + 12t^2, t\in[-1,0];
% this will compute b0 for any kernel choice provided by user
num_t1 = quad(@evalkspdf_num,-1,1,1e-14,[],kernel);
num_t2 = 0 - (-144/5 + 432/4 - 468/3 + 216/2 - 36);
denom_t1 = 0 - (-2 + 18/4 - 12/5);
denom_t2 = quad(@evalkspdf_denom,-1,1,1e-14,[],kernel);
b0 = ((num_t1^2*num_t2)/(denom_t1^2*denom_t2))^(1/5);


% [gB_bc,hb_l] = kspdf_bc(b,kernel,b,b0,usehardle);
% % right-boundary correction: reflect bids over zero and do left-boundary
% % correction 
% [gB_bc_right,hb_r] = kspdf_bc(-b,kernel,-b,b0,usehardle);
% bind = find(b >= bmax - hb_r);
% gB_bc(bind) = gB_bc_right(bind);
% GB = kscdf(b,'edf',b);
% v_bc = b' + GB./(gB_bc*(N - 1));
% v_bc = v_bc';
% vmin_bc = min(v_bc);
% vmax_bc = max(v_bc);

data = readtable('data/DataProcessed.csv');
data.PercentOfEstimates = data.PercentOfEstimates .* 100;
% data(find(data.PercentOfEstimates < 20), :) = []; 
% data(find(data.PercentOfEstimates > 200), :) = []; 



% parameters 
v = [];
bids = [];
vmin = 999;
vmax = -999;
participant_count = unique(data.ParticipantsCount);


% % step 1: recover pseudo-values using boundary correction on left boundary
for i = 1:length(participant_count)
    N = participant_count(i); % number of participants
    data_subset = data(data.ParticipantsCount == N, :);
    b = data_subset.PercentOfEstimates;
    
    if length(b) < 32
        continue 
    end 
    
    bmax = max(b);
    bmin = min(b);
    [gB_bc,hb_l] = kspdf_bc(b,kernel,b,b0,usehardle);
    % right-boundary correction: reflect bids over zero and do left-boundary
    % correction 
    [gB_bc_right,hb_r] = kspdf_bc(-b,kernel,-b,b0,usehardle);
    bind = find(b >= bmax - hb_r);
    gB_bc(bind) = gB_bc_right(bind);
    GB = kscdf(b,'edf',b);
    v_bc = b - (gB_bc)./((1-GB)*(N - 1));

    v = vertcat(v, v_bc);
    bids = vertcat(bids, b);
    
end 

v_bc = v; 


% remove inf rows
inf_idx = find(isinf(v_bc));
v_bc(inf_idx) = []; 
bids(inf_idx) = []; 
% gB_bc_global(inf_idx) = [];


vmax_bc = max(v_bc);  


% optional valuation points to evaluate each estimator at
neval = 1000;
evalpts = linspace(min(v_bc),max(v_bc),neval);

% figure 
% set(gcf,'DefaultLineLineWidth',lwidth)
% set(gca,'FontSize',fsize)
% scatter(v_bc, data.PercentOfEstimates, '.b')
% xlabel('$\mathbf{\hat{c}}$')
% ylabel('$\mathbf{b}$')
% box on

% step 2: recover valuation densities using boundary correction
[fV_bc,hv_l] = kspdf_bc(v_bc,kernel,evalpts,b0,usehardle);
% right-boundary correction: reflect pseudo-valuations over zero and do
% left-boundary correction 
[fV_bc_right,hv_r] = kspdf_bc(-v_bc,kernel,-evalpts,b0,usehardle);
vind = find(evalpts >= vmax_bc - hv_r);
fV_bc(vind) = fV_bc_right(vind);

FV = kscdf(v_bc, 'edf', evalpts);
  
if plotind == 1
    % scatter plot of bid function
    figure
    set(gcf,'DefaultLineLineWidth',lwidth)
    set(gca,'FontSize',fsize)
    scatter(v_bc, bids,'.b')
    xlabel('$\mathbf{\hat{c}}$')
	ylabel('$\mathbf{b}$')
    box on

    % bid density
%     figure
%     set(gcf,'DefaultLineLineWidth',lwidth)
%     set(gca,'FontSize',fsize)
%     plot(evalpts_b,gB_bc_global)
% 	xlabel('$\mathbf{b}$')
% 	ylabel('$\mathbf{\hat{G}_B(b)}$')
% 	box on
    
    % valuation density
    figure
    set(gcf,'DefaultLineLineWidth',lwidth)
    set(gca,'FontSize',fsize)
    plot(evalpts,fV_bc)
	xlabel('$\mathbf{\hat{c}}$')
	ylabel('$\mathbf{\hat{f}_C(c)}$')
	box on
    
    % valuation distribution
    figure 
    set(gcf, 'DefaultLineLineWidth', lwidth)
    set(gca, 'FontSize', fsize)
    plot(evalpts, FV)
    xlabel('$\mathbf{\hat{c}}$')
    ylabel('$\mathbf{\hat{F}_C(c)}$')
    box on
end


evalpts = linspace(min(v_bc), max(v_bc), max(v_bc)-min(v_bc));
FV = kscdf(v_bc, 'edf', evalpts);

if simulind == 1
    b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    prob = rand(1000, 10);
    N = 16;
    
    
    
    for i = 1:1000
        for j = 1:10 
            idx = find(prob(i, j) == FV);
            
            if isempty(idx) 
                [minValue, closestIdx] = min(abs(FV - prob(i, j))); 
                cost = evalpts(closestIdx);
                fprintf('%d %d\n', i, j)
                b(j) = b(j) + cost + ... 
                    trapz(evalpts(closestIdx:end), 1-FV(closestIdx:end)) / (1-FV(closestIdx))^N;
            else 
                cost = mean(evalpts(idx));
                b(j) = b(j) + cost + ... 
                    trapz(evalpts(idx(1):end), 1-FV(idx(1):end)) / (1-prob(i, j))^N;
            end 
        end 
    end 
    
    b = b ./ 1000;
end 