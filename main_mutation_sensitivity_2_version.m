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
    
    MIM_mut=f_compute_eq_mutated_CRN_new_way(protein(i), new_CMIM, idx_basic_species, x_0_phys, k_values);
    
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
    ris_mut.(protein{i}).react_2_be_removed=MIM_mut.info.react_rem_lof;
    ris_mut.(protein{i}).idx_basic_species=MIM_mut.species.idx_basic_species;
    rho=MIM_mut.matrix.rho;
    rho(MIM_mut.species.idx_law_rem)=0;
    % 
     [ris.(protein{i}).SSI_k, ris.(protein{i}).SSI_c]=f_compute_SSI(idx_sp, x_eq, MIM_mut.rates.std_values,...
         ris_mut.(protein{i}).S, MIM_mut.matrix.Nl, rho, ...
         MIM_mut.species.idx_basic_species, MIM_mut.matrix.v);
    clear MIM_mut


    [ris.(protein{i}).SSI_k_sort, ris.(protein{i}).ord_SSI_k]=sort(ris.(protein{i}).SSI_k, 'descend');

end