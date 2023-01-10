function [x0_mut, rateconst_mut] = f_define_mutated_condition(protein, ...
                                    x0, rateconst, MIM, perc)

%% Input:
% protein : string
%   Protein affected by mutation. Only a limited number of protein are
%   allowed (see variable all_proteins)
% x0 : float [n_species x 1]
%   Initial condition before mutation
% rateconst : float [n_reactions x 1]
%   Rate constants before mutation
% MIM : structure
%   The MIM
% perc : float
%   degree of LoF mutation. (perc = 0 --> null mutation, i.e. the concentration 
%   of the mutated protein is set to 0)
%   By now this variable should be let empty in GoF


%% Output
% x0_mut : float [n_species x 1]
%  Initial condition after mutation
% rateconst_mut : float [n_reactions x 1]
%  Rate consants after mutation.

%% Step 1. Define implemented protein
lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_proteins = [lof_mutations, gof_mutations, lof_mutations_type2];

if ~ismember(protein, all_proteins)
    error('Mutation for %s has not been implemented yet', protein)
end

if nargin == 4
    perc = 0;
end

%% Step 2. Define general variable
[~, idx_protein] = ismember(protein, MIM.species.names);

Nl = MIM.matrix.Nl;
cons_law_protein = find(Nl(:, idx_protein));
species_cons_law_protein = find(Nl(cons_law_protein, :));

idx_basic_species = find(MIM.species.std_initial_values~=0);

if isempty(cons_law_protein)
   fprintf('%s does not belong to any conservation law \n', protein)
else
   fprintf('Species involved in %s''s conservation law: \n', protein)
   disp(MIM.species.names(species_cons_law_protein))
end


%% Step 3. Define parameters in the mutated cell
% NOTE: A loss of function (LoF) mutation corresponds to a mutated initial
% condition, while a gain of function (GoF) mutation corresponds to 
% different values of the constant rates.

%% 3.1. Initial condition
switch protein
    
    case lof_mutations % Loss of function mutation
        
        % 3.1.1. Project the values of the total concentrations
        rho = Nl*x0;
        rho(cons_law_protein) = perc*rho(cons_law_protein);

        % 3.1.2. Decompose Nl
        Nl_1 = Nl(:, idx_basic_species);
        Nl_2 = Nl; Nl_2(:, idx_basic_species) = [];

        % 3.1.3. Intialize everything to zero
        x0_mut = zeros(size(x0));

        % 3.1.4 Species not element of the network and not in CL
        idx_nobasic_species = 1:numel(x0); 
        idx_nobasic_species(idx_basic_species) = [];
        aux_sp = setdiff(idx_nobasic_species, species_cons_law_protein);
        x0_mut(aux_sp) = x0(aux_sp);

        % 3.1.5 Elements of the network
        x0_mut(idx_basic_species) = ...
            Nl_1 \ (rho - Nl_2 * x0_mut(idx_nobasic_species));

    case gof_mutations % Gain of function mutations
        
        x0_mut = x0;
    
    case lof_mutations_type2

        switch protein

            case 'TP53' % TP53

                dead_sps = 'TP53_generator'; 
        end

        [~, idx_dead_sps] = ismember(dead_sps, MIM.species.names);
        x0_mut = x0; x0_mut(idx_dead_sps) = 0;

end

%% 2.b. Rate constants
switch protein
    
    case [lof_mutations, lof_mutations_type2] % Loss of function mutation or TP53
        
    rateconst_mut = rateconst;    
        
    case gof_mutations % Gain of function mutations
    
    switch protein
    
        case 'Raf'
            idx_reactions = [68, 70];
        case 'Ras'
            idx_reactions = [97, 99, ...
                        125, 127, ...
                        128, 130, ...
                        321, 323, ...
                        346, 348, ...
                        349, 351, ...
                        448, 450, ...
                        473, 475, ...
                        476, 478, ...
                        44];
       case 'PI3K'
            idx_reactions = [181, 183, ...
                             184, 186, ...
                             382, 384, ...
                             385, 387, ...
                             505, 507, ...
                             508, 510, ...
                             593, 595, ...
                             596, 598, ...
                             599, 601, ...
                             686, 688, ...
                             178, 180];
        case 'BetaCatenin'
            idx_reactions = [717, 719, ...
                             773, 775];
    end
    
    % Some sanity checks
    original_S = MIM.matrix.S;
    mutated_S = MIM.matrix.S; mutated_S(:, idx_reactions) = 0;
    if rank(mutated_S) ~= rank(original_S)
        error('%s does not satisfy the condition for gain of function ', protein)
    end
    
    rateconst_mut = rateconst;
    rateconst_mut(idx_reactions) = perc * rateconst_mut(idx_reactions);
        
end

end
