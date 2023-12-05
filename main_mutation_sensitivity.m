%% Clear all
clc
clear
close all

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

addpath(fullfile('.', 'funcs'))
folder_results = fullfile('.', 'results');
if ~exist(folder_results, 'dir')
    mkdir(folder_results)
end
warning('off', 'all')

%% Step 1. Load data
% Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug.mat'); % Da cambiare in base a dove metti i file
load(file_crn, 'new_CMIM')
file_crc_crn = fullfile('data', 'CRC_CRN.mat');
load(file_crc_crn, 'CMIM')
list_reactions = CMIM.reactions; clear CRC

%% Step 2. Import the initial data
% number of species, reactions, conservation laws
n_species = numel(new_CMIM.species.names);
n_reactions = size(new_CMIM.matrix.S, 2);
n_cons_laws = size(new_CMIM.matrix.Nl, 1);


% define the new matrices and vectors of the system
Sm = new_CMIM.matrix.S;
vm = new_CMIM.matrix.v;
k_values=new_CMIM.rates.std_values;
cons_laws=new_CMIM.matrix.Nl;
x_0_phys=new_CMIM.species.std_initial_values;
idx_basic_species = find(x_0_phys>0);
max_counter=300;
idx_one=new_CMIM.matrix.ind_one;
rho_phys=cons_laws*x_0_phys;
ris_phys=f_NLPC_restart(x_0_phys, k_values, Sm, cons_laws, rho_phys,...
    idx_basic_species, vm,idx_one, max_counter, 0);
x_eq_phys = ris_phys.x;

x_values=x_eq_phys;
idx_sp_phys=1:numel(x_values);
rho_phys=cons_laws*x_0_phys;
[SSI_k_phys, SSI_c_phys]=f_compute_SSI(idx_sp_phys, x_values, k_values, ...
    Sm, cons_laws, rho_phys, idx_basic_species, vm);


%% Step 3. Create the mutatated CRC
lof_mutation = {'APC', 'SMAD4'};
gof_mutation = {'Ras'};
lof_mutation_type2 = {'TP53'};
all_proteins = [gof_mutation, lof_mutation, lof_mutation_type2];

protein=["Ras", "APC", "SMAD4", "TP53"];

for i=1:numel(protein)
    
    MIM_mut=f_compute_eq_mutated_CRN(protein(i), new_CMIM, idx_basic_species, x_0_phys, k_values);
    
    ris_mut.(protein{i}).x_eq=MIM_mut.species.x_eq;
    null_species=MIM_mut.species.null_species;
    sp_rem=MIM_mut.species.sp_rem;
    sp_to_rem=sort([null_species, sp_rem]);
    ris_mut.(protein{i}).react_rem=MIM_mut.info.react_rem;
    %% Sensitivity matrices for the mutated CRC-CRN
    x_eq=ris_mut.(protein{i}).x_eq;
    x_eq(null_species)=0;
    idx_sp=find(x_eq);
    ris_mut.(protein{i}).n_basic_mut=MIM_mut.species.n_basic_mut;
    ris_mut.(protein{i}).S=MIM_mut.matrix.S;
    ris_mut.(protein{i}).idx_law_rem=MIM_mut.species.idx_law_rem;
    
    ris_mut.(protein{i}).idx_basic_species=MIM_mut.species.idx_basic_species;

    [ris.(protein{i}).SSI_k, ris.(protein{i}).SSI_c]=f_compute_SSI(idx_sp, x_eq, MIM_mut.rates.std_values,...
        ris_mut.(protein{i}).S, MIM_mut.matrix.Nl, MIM_mut.matrix.rho, ...
        MIM_mut.species.idx_basic_species, MIM_mut.matrix.v);
    clear MIM_mut
end
%% Graphs
% SSI on k
f_eff_mut_log_k = figure('units','normalized','outerposition',[0 0  0.75 1]);

for i=1:numel(protein)
    aux_SSI_k_phys=SSI_k_phys;
    aux_SSI_k=ris.(protein{i}).SSI_k;
    react_rem=ris_mut.(protein{i}).react_rem;
    % delete less sensible reactions for KRAS
    if string(protein(i))=="Ras"
        disp(protein(i))
        react_rem=ris_mut.Ras.react_rem;
        ind_react=1:1:n_reactions;
        ind_react_mut=ind_react;
        ind_react_mut(react_rem)=[];
        less_sens=find(ris.Ras.SSI_k<10^-30);
        react_less=numel(less_sens);
        ind_orig_less_sens=ind_react_mut(less_sens);
        disp(list_reactions.arrow(list_reactions.reactions2flux_rates(less_sens)))
        react_rem=sort([react_rem;ind_orig_less_sens']);
        aux_SSI_k(less_sens)=[];
        disp(react_rem)
    end
    disp(react_rem)
    aux_SSI_k_phys(react_rem)=[];
    [aux_SSI_k_phys_sort, order_phys]=sort(aux_SSI_k_phys, 'descend');
    
    
    subplot(2,2,i)
    
    semilogy(aux_SSI_k(order_phys), 'k', 'linewidth', 3)
    hold on
    semilogy(aux_SSI_k_phys_sort, 'r--', 'linewidth', 3)
    ylim([10^-16, 10^0])
    set(gca, 'Fontsize', 15)
    
    switch protein(i)
        case "Ras"
            aux_title = 'GoF KRAS (a)';
        case "APC"
            aux_title = 'LoF APC (b)';
        case "SMAD4"
            aux_title = 'LoF SMAD4 (c)';
        case "TP53"
            aux_title = 'LoF TP53 (d)';
    end
    
    switch i
        case 1
            ylabel('e^k_j')
            
        case 3
            ylabel('e^k_j')
            xlabel('index j')
            
        case 4
            xlabel('index j')
            
    end
    
    
    xlim([1 size(ris_mut.(protein{i}).S,2)-numel(react_rem)])
    xticks(0:200:size(ris_mut.(protein{i}).S,2)-numel(react_rem)) % controllare
    title(aux_title)
    
    lg = legend('Mutated', 'Phys', 'Location', 'southwest');
    lg.FontSize = 18;
    
end

% SSI on c
f_eff_mut_log_c = figure('units','normalized','outerposition',[0 0  0.75 1]);

for i=1:numel(protein)
    aux_SSI_c_phys=SSI_c_phys;
    aux_SSI_c=ris.(protein{i}).SSI_c;
    aux_SSI_c_phys(ris_mut.(protein{i}).idx_law_rem)=[];
    [aux_SSI_c_phys_sort, order_phys]=sort(aux_SSI_c_phys, 'descend');
    
    
    subplot(2,2,i)
    
    semilogy(aux_SSI_c(order_phys), 'k', 'linewidth', 3)
    hold on
    semilogy(aux_SSI_c_phys_sort, 'r--', 'linewidth', 3)
    set(gca, 'Fontsize', 15)
    
    switch protein(i)
        case "Ras"
            aux_title = 'GoF KRAS (a)';
        case "APC"
            aux_title = 'LoF APC (b)';
        case "SMAD4"
            aux_title = 'LoF SMAD4 (c)';
        case "TP53"
            aux_title = 'LoF TP53 (d)';
    end
    
    switch i
        case 1
            ylabel('e^c_j')
        case 3
            ylabel('e^c_j')
            xlabel('index j')
        case 4
            xlabel('index j')
    end
    
    title(aux_title)
    xlim([1 ris_mut.(protein{i}).n_basic_mut])
    xticks(0:20:ris_mut.(protein{i}).n_basic_mut)
    ylim([10^-4, 1.1*10^-1])
    
    lg = legend('Mutated', 'Phys', 'Location', 'northeast');
    lg.FontSize = 18;
   
end
saveas(f_eff_mut_log_c, fullfile(folder_results, 'single_gene_mutations_c.png'))
%% Table for KRAS
SSI_k_phys(ris_mut.Ras.react_rem)=[];
[values_abs_rel_diff, order_diff]=sort(abs((ris.Ras.SSI_k- SSI_k_phys)./SSI_k_phys), 'descend');
values_rel_diff=(ris.Ras.SSI_k-SSI_k_phys)./SSI_k_phys;
ind_react(ris.Ras.react_rem)=[]; %tolgo le reazioni rimosse

table(SSI_k_phys(order_diff), ris.Ras.SSI_k(order_diff), values_rel_diff(order_diff))
table(string(list_reactions.details(ind_react(order_diff),1)), values_rel_diff(order_diff), 'VariableNames', {'Reactions', 'RelDiffSSI'})

table_file=fullfile(folder_results, 'reactions_Ras.txt');
fileID = fopen(table_file, 'w');
disp(['Writing on ', table_file, '...'])

% row_table=10;
% for ii = 1:row_table
%     fprintf(fileID, '%s & \\verb| %s | & %1.2e  \\\\ \\hline \n', ...
%         strcat('R', num2str(order_diff(ii))), ...
%         string(list_reactions.details(order_diff(ii),1)), ...
%         (values_rel_diff((order_diff(ii)))));
% end
% fclose(fileID);



row_table=10;
for ii = 1:row_table
    fprintf(fileID, '%s & \\verb| %s | & %1.2e  \\\\ \\hline \n', ...
        strcat('R', num2str(order_diff(ii))), ...
        string(list_reactions.details(order_diff(ii),1)), ...
        (values_rel_diff));
end
fclose(fileID);
