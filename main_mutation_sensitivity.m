%% Clear all

%% unica threshold -> uniche specie che tolgo sono quelle delle gof e delle lof

clc
clear
close all

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

addpath(fullfile('.', 'funcs'))
folder_results = '.\results';
if ~exist(folder_results, 'dir')
    mkdir(folder_results)
end
warning('off', 'all')

%% 1. Load data
% Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug_complete.mat'); % Da cambiare in base a dove metti i file
load(file_crn, 'CMIM')

%% 2. Import the initial data
% number of species, reactions, conservation laws
n_species = numel(CMIM.species.names);
n_reactions = size(CMIM.matrix.S, 2);
n_cons_laws = size(CMIM.matrix.Nl, 1);

% define the new matrices and vectors of the system
list_reactions = CMIM.reactions;
Sm = CMIM.matrix.S;
vm = CMIM.matrix.v;
k_values=CMIM.rates.std_values;
cons_laws=CMIM.matrix.Nl;
x_0_phys=CMIM.species.std_initial_values;
idx_basic_species = find(x_0_phys>0);
max_counter=300;
idx_one=CMIM.matrix.ind_one;
rho_phys=cons_laws*x_0_phys;
ris_phys=f_NLPC_restart(x_0_phys, k_values, Sm, cons_laws, rho_phys,...
    idx_basic_species, vm,idx_one, max_counter, 0);
x_values = ris_phys.x;

toll=10^-15;
%% 3. Physiological status

idx_sp_phys=1:numel(x_values);
rho_phys=cons_laws*x_0_phys;
[SSI_k_phys, SSI_c_phys]=f_compute_SSI(idx_sp_phys, x_values, k_values, ...
    Sm, cons_laws, rho_phys, idx_basic_species, vm);


%% 4. Create the mutatated CRC
lof_mutation = {'APC', 'SMAD4'};
gof_mutation = {'Ras'};
lof_mutation_type2 = {'TP53'};
all_proteins = [gof_mutation, lof_mutation, lof_mutation_type2];

protein=["Ras", "APC", "SMAD4", "TP53"];


f_thresh = figure('units','normalized','outerposition',[0 0  0.75 1]);
for i=1:numel(protein)

    MIM_mut=f_compute_eq_mutated_CRN(protein(i), CMIM, idx_basic_species, x_0_phys, k_values);

    ris_mut.(protein{i}).x_0=MIM_mut.species.x_0;
    ris_mut.(protein{i}).x_eq=MIM_mut.species.x_eq;
    ris_mut.(protein{i}).react_rem=MIM_mut.info.react_rem;
    ris_mut.(protein{i}).S=MIM_mut.matrix.S;
    ris_mut.(protein{i}).idx_law_rem=MIM_mut.species.idx_law_rem;
    ris_mut.(protein{i}).idx_basic_species=MIM_mut.species.idx_basic_species;
    ris_mut.(protein{i}).null_species=MIM_mut.info.null_species;
    ris_mut.(protein{i}).n_basic_mut=MIM_mut.info.n_basic_mut;
    ris_mut.(protein{i}).react_rem_lof=MIM_mut.info.react_rem_lof;
    x_eq=ris_mut.(protein{i}).x_eq;
    
    %x_eq(ris_mut.(protein{i}).null_species)=0;
    idx_sp=find(x_eq>toll);
    idx_sp_rem=setdiff(1:n_species, idx_sp);
    idx_sp_no_zero=find(ris_mut.(protein{i}).x_eq);
    
    subplot(2,2,i)
    semilogy(sort(ris_mut.(protein{i}).x_eq(idx_sp_no_zero), 'descend'), '-o')
    xlim([1, numel(idx_sp_no_zero)])
    xlabel('species')
    ylabel('concentration at equilibrium')

    disp("*Removing species...*")
    disp(CMIM.species.names(idx_sp_rem))

    [ris.(protein{i}).SSI_k, ris.(protein{i}).SSI_c]=f_compute_SSI(idx_sp, x_eq, MIM_mut.rates.std_values,...
        ris_mut.(protein{i}).S, MIM_mut.matrix.Nl, MIM_mut.matrix.rho, ...
        MIM_mut.species.idx_basic_species, MIM_mut.matrix.v);

    % clear MIM_mut
end
saveas(f_thresh, fullfile(folder_results, 'choose_threshold.png'))
%% Graphs
% SSI on k
f_eff_mut_log_k = figure('units','normalized','outerposition',[0 0  0.75 1]);

for i=1:numel(protein)
    aux_SSI_k_phys=SSI_k_phys;
    aux_SSI_k=ris.(protein{i}).SSI_k;
    react_rem=ris_mut.(protein{i}).react_rem;
    react_rem_2=ris_mut.(protein{i}).react_rem_lof;
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
    end
    aux_SSI_k_phys(react_rem)=[];
    aux_SSI_k_phys(react_rem_2)=[];
    aux_SSI_k(react_rem_2)=[];
    [aux_SSI_k_phys_sort, order_phys]=sort(aux_SSI_k_phys, 'descend');


    subplot(2,2,i)

    semilogy(aux_SSI_k(order_phys), 'k', 'linewidth', 3)
    hold on
    semilogy(aux_SSI_k_phys_sort, 'r--', 'linewidth', 3)
    ylim([10^-30, 10^0])
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
            ylabel('e^k_j','FontSize', 20, 'FontName', 'Italic')

        case 3
            ylabel('e^k_j','FontSize', 20, 'FontName', 'Italic')
            xlabel('index j','FontSize', 20, 'FontName', 'Italic')

        case 4
            xlabel('index j','FontSize', 20, 'FontName', 'Italic')

    end

    xlim([1 size(ris_mut.(protein{i}).S,2)-max(numel(react_rem), numel(react_rem_2))])
    xticks(0:200:size(ris_mut.(protein{i}).S,2)-max(numel(react_rem), numel(react_rem_2))) % controllare
    title(aux_title,'FontSize', 20)

    lg = legend('Mutated', 'Phys', 'Location', 'southwest');
    lg.FontSize = 20;

end

saveas(f_eff_mut_log_k, fullfile(folder_results, 'single_gene_mutations_k.png'))
% SSI on c
f_eff_mut_log_c = figure('units','normalized','outerposition',[0 0  0.75 1]);

for i=1:numel(protein)
    aux_SSI_c_phys=SSI_c_phys;
    aux_SSI_c=ris.(protein{i}).SSI_c;
    aux_SSI_c(ris_mut.(protein{i}).idx_law_rem)=[];
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
            ylabel('e^c_j','FontSize', 20, 'FontName', 'Italic')
        case 3
            ylabel('e^c_j','FontSize', 20, 'FontName', 'Italic')
            xlabel('index j', 'FontSize', 20, 'FontName', 'Italic')
        case 4
            xlabel('index j','FontSize', 20, 'FontName', 'Italic')
    end

    title(aux_title, 'FontSize', 20, 'FontName', 'Italic')
    xlim([1 ris_mut.(protein{i}).n_basic_mut])
    xticks(0:20:ris_mut.(protein{i}).n_basic_mut)
    ylim([10^-4, 1.1*10^-1])

    lg = legend('Mutated', 'Phys', 'Location', 'northeast');
    lg.FontSize = 20;
   

end
saveas(f_eff_mut_log_c, fullfile(folder_results, 'single_gene_mutations_c.png'))
%% Table for KRAS
ind_react=1:1:n_reactions;
SSI_k_phys_saved=SSI_k_phys;
SSI_k_phys(ris_mut.Ras.react_rem)=[];
[values_abs_rel_diff, order_diff]=sort(abs((ris.Ras.SSI_k- SSI_k_phys)./SSI_k_phys), 'descend');
values_rel_diff=(ris.Ras.SSI_k-SSI_k_phys)./SSI_k_phys;
ind_react(ris_mut.Ras.react_rem)=[]; %tolgo le reazioni rimosse

table(SSI_k_phys(order_diff), ris.Ras.SSI_k(order_diff), values_rel_diff(order_diff))
table(string(list_reactions.details(ind_react(order_diff),1)), values_rel_diff(order_diff), 'VariableNames', {'Reactions', 'RelDiffSSI'})

table_file=fullfile(folder_results, 'reactions_Ras.txt');
fileID = fopen(table_file, 'w');
disp(['Writing on ', table_file, '...'])

row_table=10;
for ii = 1:row_table
    fprintf(fileID, '%s & \\verb| %s | & %1.2e  \\\\ \\hline \n', ...
        strcat('R', num2str(order_diff(ii))), ...
        string(list_reactions.details(order_diff(ii),1)), ...
        (values_rel_diff((order_diff(ii)))));
end
fclose(fileID);

