clc
clear
close all

rand('seed', abs(round((now * 100000 - round(now*100000))*1000)));

folder_code = fullfile('.', 'func');
addpath(folder_code)

%% Step 1. Define input data
% 1.1. Folders
target_folder = fullfile('.', 'data');
folder_results = fullfile('.', 'results');

% 1.2. MIM matrices
file_mim = fullfile(target_folder, 'AML_CRN.mat');

aux_file_ris ='ris_uniqueness_rho_%d.mat';

% 1.4. Load 
load(file_mim, 'CMIM');

% 1.5. Parameters for the simulation
%   - Number of runs
n_rho = 100;
n_sim = 30; % --> ****** Increase this number when using Newton *******

%   - Parameters of CMIM model
rate_constants = CMIM.rates.std_values;
max_t = 10^8;

%   - auxiliar variables for the simulation
true_x0 = CMIM.species.std_initial_values;
idx_basic_species = find(true_x0 >0);
Nl = CMIM.matrix.Nl;

[n_cl, n_species] = size(Nl);
idx_no_cl = find(sum(Nl, 1)==0);
idx_in_cl = 1:n_species; idx_in_cl(idx_no_cl) = [];
idx_sp_todo = setdiff(idx_in_cl, idx_basic_species);

bound_sp_nocl = 10^4;
bound_sp = 1;

par_loguni_rho = [-2, 3];
par_loguni_sp_nocl = [-5, 5];

%% Corr cl
% corr_cl = zeros(n_cl, n_cl);
% for ir = 1:n_cl
%     for jr = 1:n_cl
%         corr_cl(ir, jr) = numel(intersect(find(Nl(ir, :) > 0), find(Nl(jr, :)>0)));
%     end
%      corr_cl(ir, ir) = 0;
% end
% 
% imagesc(corr_cl)
% colorbar

%% Step 2. Run simulation

%%   2.1. Initialize output variables;
all_x0 = zeros(n_species, n_sim);
all_xe = zeros(n_species, n_sim);
all_ntimes = zeros(n_sim, 1);
all_dx_eq = zeros(n_sim, 1);

aux_all_rho = 10.^((par_loguni_rho(2) - par_loguni_rho(1))*rand(n_cl, n_rho) ...
                +par_loguni_rho(1));

for rr = 1%:n_rho
    
    %% 2.2. Draw the stoichiometric surface
    if rr  == 1
        rho = Nl*true_x0;
    else
        rho = aux_all_rho(:, rr);
    end
    
    file_ris = fullfile(folder_results, sprintf(aux_file_ris, rr));
     
for ii = 1:n_sim 

    fprintf('Running simulation num = %d for sup num = %d \n', ii, rr)
    
%%  2.2. Draw initial value
    x0 = zeros(n_species, 1);
    
    % - Species that do not belong to any conservation law
    x0(idx_no_cl) = 10.^((par_loguni_sp_nocl(2) - par_loguni_sp_nocl(1))...
        *rand(numel(idx_no_cl), 1)+par_loguni_sp_nocl(1));
    
    % - Species in cls but not basic
    ordered_sp = idx_sp_todo(randperm(length(idx_sp_todo)));
    aux_rho = rho;
    for is = ordered_sp
        cl_is = find(Nl(:, is)>0);
        tmp_value = bound_sp * rand(1) * min(aux_rho(cl_is) ./ Nl(cl_is, is));
        x0(is) = tmp_value;
        aux_rho(cl_is) = aux_rho(cl_is) - Nl(cl_is, is) * tmp_value;
        if any(aux_rho<0)
            disp(is)
        end
    end
    
    %% SARA
    x0(67) = 0; x0(68) = 0;
    
    % - Basic species
    for ir = 1:n_cl
        tmp_basic = intersect(idx_basic_species, find(Nl(ir, :)));
        x0(tmp_basic) = max(rho(ir) - Nl(ir, :)*x0, 0);
    end
    
    fprintf('||C_t - C_t_vero|| = %2.4e \n', norm(rho - Nl*x0, 'Inf'))
    fprintf('min(x_0) = %2.2e \n', min(x0))
    
%%  2.3. Find corresponding equilibrium
    [time, x_t] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rate_constants, CMIM, 'Sv'), [0 max_t], x0);
    x_t = x_t'; % Size: n_species x n_times
    n_times = numel(time); 
    xe = x_t(:, end);
    
    fprintf('||C_t_rec - C_t_vero|| = %2.4e \n', norm(rho - Nl*xe, 'Inf'))
    
%%  2.4 Store and clear variables
    all_x0(:, ii) = x0;
    all_xe(:, ii) = xe;
    all_ntimes(ii) = n_times;
    all_dx_eq(ii) = norm(f_odefun_MIM(0, xe, rate_constants, CMIM, 'Sv'), 'inf');
    
%     figure
%     plot(time, x_t(67, :), 'linewidth', 2)
%     pause
    
    clear aux_elements x0 time x_t n_time xe
    
end

%% Step 3. Save results
disp('Saving results...')
ris_uni.all_x0 = all_x0;
ris_uni.all_xe = all_xe;
ris_uni.all_ntimes = all_ntimes;
ris_uni.all_dx_eq = all_dx_eq;
ris_uni.rho = rho;
save(file_ris, 'ris_uni')

end



