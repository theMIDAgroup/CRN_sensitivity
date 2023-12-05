clc
clear
close all

% This code creates the figures representing the correctness of the SSI

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 
warning('off', 'all')

addpath('./funcs')

%% Step 1. Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug_complete.mat');
load(file_crn)

%% Step 2. Import the initial data
vm = CMIM.matrix.v;
n_species = numel(CMIM.species.names);
n_reactions = size(CMIM.matrix.S, 2);
n_cons_laws = size(CMIM.matrix.Nl, 1);
ind_one = CMIM.matrix.ind_one;

k_values = CMIM.rates.std_values;
Sm = CMIM.matrix.S;
x_0 = CMIM.species.std_initial_values;
idx_basic_species = find(x_0>0);
cons_laws = CMIM.matrix.Nl;
rho=cons_laws*x_0;

mut_lof=0;
%% Step 3. Compute the values of the SSI

%% 3.1. Compute equilibrium
max_counter = 300;
%ris_phys = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
%           vm, ind_one, max_counter, 0);
% x_values = ris_phys.x;
file_x_phys = fullfile('results', 'x_phys.mat');
load(file_x_phys)
% save("./results/x_phys.mat","x_values")
%% 3.2. Analytically compute the sensitivity matrix
% Compute the analytic jacobian of v 
idx_sp=1:numel(x_values);
[SSI_k, SSI_c]=f_compute_SSI(idx_sp, x_values, k_values, ...
                                 Sm, cons_laws, rho, idx_basic_species, vm);
%% 4. Tables
% Table SSI_k

[SSI_k_sorted, index_SSI_k_sorted]=sort(SSI_k, 'descend');
arrows_reactions = CMIM.reactions.details(:,1);
T_SSI_k=table(SSI_k_sorted, string(arrows_reactions(index_SSI_k_sorted)), 'VariableNames', {'e^k_j', 'j-th reaction'});

table_file=fullfile('results', 'reactions_ordered_by_SSI.txt');
fileID = fopen(table_file, 'w');
disp(['Writing on ', table_file, '...'])

row_table=numel(SSI_k);
for ii = 1:row_table
    fprintf(fileID, '%s & \\verb| %s | & %1.2e  \\\\ \\hline \n', ...
        strcat('R', num2str(index_SSI_k_sorted(ii))), ...
        string(arrows_reactions(index_SSI_k_sorted(ii))), ...
        SSI_k_sorted(ii));
end
fclose(fileID);



% Table SSI_k

[SSI_c_sorted, index_SSI_c_sorted]=sort(SSI_c, 'descend');

colored_species_alias = CMIM.species.alias_conc;
table_cons_laws=fullfile('results', 'conservation_laws_ordered_by_SSI.txt');
fileID = fopen(table_cons_laws, 'w');
disp(['Writing on ', table_cons_laws, '...'])
for ii = 1:size(cons_laws, 1)
    i_s=index_SSI_c_sorted(ii);
    idx_sp = find(cons_laws(i_s,:));
    % string = '$';
    % for jj = 1:numel(idx_sp)
    %     string = strcat(string , ...
    %         ['+ ', num2str(cons_laws(i_s, idx_sp(jj))), ...
    %         '\\ ', colored_species_alias{idx_sp(jj)}]);
    % end
    % string = strcat(string, ' $ = $\rho_{', num2str(i_s), '} $ & | %1.2e ', num2str(SSI_c_sorted(ii)), '\\\\ \hline \n');
    % fprintf(fileID, string);
    string='';
    for jj = 1:numel(idx_sp)
        disp(jj)
        string = strcat(string, '+ ', num2str(cons_laws(i_s, idx_sp(jj))), ...
            '\\ ', colored_species_alias{idx_sp(jj)});
        disp(string)
    end
    disp(string)
    string = strcat(string , '=\rho_{', num2str(i_s), '}');
    fprintf(fileID, '$ %s $ \\ & %1.2e  \\\\ \\hline \n', ...
        string, SSI_c_sorted(ii));

end
fclose(fileID);
