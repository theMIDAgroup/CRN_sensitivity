function [CRN, N] = f_add_drug_Raf_from_file(CRN, drug)
    
    %modified fullfile('../data', sprintf('%s', drug))
    protein_folder = fullfile('data', sprintf('%s', drug));
    % Files txt with data
    drug_file_species1 = fullfile(protein_folder, 'add_drug_species_1.txt');
    drug_file_species2 = fullfile(protein_folder, 'add_drug_species_2.txt');
    drug_file_reactions1 = fullfile(protein_folder, 'add_drug_reactions_1.txt');
    
	%% Store data
    % tab delimiter
    delimiter = '\t';    
    % read columns of data
    formatSpec = '%s%f%s%f';
    formatLabel = '%s%s%s%s';
    % open the text file.
    fileID = fopen(drug_file_species1,'r');
    % read columns of data according to format string.
    labels_sp1 = textscan(fileID, formatLabel, 1, 'Delimiter', delimiter,  'ReturnOnError', false);
    file_species1 = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false, 'headerlines', 1);
    % close the text file
    fclose(fileID); 
    
    formatSpec = '%s%s';
    fileID = fopen(drug_file_species2,'r');
    file_species2 = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false, 'headerlines', 1);
    fclose(fileID);
    
    formatSpec = '%s%s%s';
    fileID = fopen(drug_file_reactions1,'r');
    file_reactions1 = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false, 'headerlines', 1);
    fclose(fileID);
    
    %% Add drugs
    % 3.a. Species
    ind_names = find(strcmp([labels_sp1{:}], 'names'));
    ind_inval = find(strcmp([labels_sp1{:}], 'std_in_value'));
    ind_unit = find(strcmp([labels_sp1{:}], 'units'));
    ind_const = find(strcmp([labels_sp1{:}], 'is_constant'));
    n_sp = numel(CRN.species.names);    
    n_prod = length(CRN.species.products);
    N = length(file_species1{1});
    for i=1:N
        if ismember(char(file_species1{ind_names}(i)), CRN.species.names)
            error('This specie already exists');
        end
        CRN.species.names{n_sp+i} = char(file_species1{ind_names}(i));
        CRN.species.std_initial_values(n_sp+i) = file_species1{ind_inval}(i);
        CRN.species.units{n_sp+i} = char(file_species1{ind_unit}(i));
        CRN.species.is_constant(n_sp+i) = file_species1{ind_const}(i); 
        % Aliases
        CRN.species.alias{n_sp+i} = strcat(['A_{', num2str(n_sp+i), '}']);
        CRN.species.alias_conc{n_sp+i} = strcat(['x_{', num2str(n_sp+i), '}']);
    end
    for j = 1:length(file_species2{1})
        for i=1:2
            % Products
            CRN.species.products{n_prod+j,i} = char(file_species2{i}(j));
            [~, CRN.species.products_ind(n_prod+j,i)] = ismember(CRN.species.products{n_prod+j,i}, CRN.species.names);
        end
    end

    % 3.b Reaction and rate constant
    aux_n_reac = numel(CRN.reactions.arrow);
    aux_n_par = numel(CRN.rates.names);
    j = 1;
    if (drug == 'DBF')
        s = j;
    elseif (drug == 'TMT')
         %s = 3;
        s = length(CRN.reactions.details) - 851 + 1;
    end
    
    for i=1:length(file_reactions1{1})
        CRN.reactions.arrow{aux_n_reac+i} = char(file_reactions1{1}(i));
        CRN.rates.in_reactions{aux_n_par + j} = CRN.reactions.arrow{aux_n_reac+i};
        CRN.rates.names{aux_n_par + j} = strcat(['cd_', num2str(s)]);
        CRN.rates.std_values(aux_n_par + j) = 0; %capire perchè 0!!
        CRN.rates.units{aux_n_par + j} = '1/(nanomole*second)';
        CRN.reactions.Flux_rate{aux_n_reac+i} = ....
        strcat([CRN.rates.names{aux_n_par+j}, char(file_reactions1{2}(i))]); %i
        j = j+1;
        s = s+1;
        if contains(char(file_reactions1{1}(i)), '<->')
            CRN.rates.names{aux_n_par + j} = strcat(['cd_', num2str(s)]);
            CRN.rates.in_reactions{aux_n_par + j} = CRN.reactions.arrow{aux_n_reac+i};
            CRN.rates.std_values(aux_n_par + j) = 0; %capire perchè 0!!
            CRN.rates.units{aux_n_par + j} = '1/second'; % qua dovrebbe andarci un +2 - SILVIA
            CRN.reactions.Flux_rate{aux_n_reac+i} = strcat([CRN.reactions.Flux_rate{aux_n_reac+i}, CRN.rates.names{aux_n_par+j},  char(file_reactions1{3}(i))]);
            j = j+1;
            s = s+1;
        end
    end
    
    %% Step 4. Decode reactions
    % 4.1. Rate constants
    [~, CRN.rates.in_reactions_idx] = ismember(CRN.rates.in_reactions, CRN.reactions.arrow);

    % 4.2. Reactions
    [CRN.reactions.Flux_decoded, CRN.reactions.Flux_decoded_ind,...
        CRN.species.products, CRN.species.products_ind,...
        complexes, complexes_ind, CRN.reactions.details, CRN.reactions.details_ind,...
        CRN.reactions.nr_irreversible, CRN.reactions.nr_reversible, ...
        CRN.reactions.reactions2flux_rates] = ...
        decode_reactions(CRN.reactions.arrow, CRN.reactions.Flux_rate, ...
                         CRN.species.names, CRN.rates.names);

    % 4.3. Rename reactions with a more mathematical notation
    [n_reactions, aux_n_cols] = size(CRN.reactions.details);
    for ir = 1:n_reactions
        CRN.reactions.details(ir, aux_n_cols+1) = {strcat(['R_{', num2str(ir),'}'])};
        CRN.rates.alias(CRN.reactions.details_ind(ir, 1)) =  ...
                            {strcat(['k_{', num2str(ir), '}'])};
    end
    
    %% Step 5. Define system in matrix form
    [CRN.matrix.S, CRN.matrix.v, CRN.matrix.Z, CRN.matrix.B, CRN.matrix.ind_one] = matrix_ZBv(CRN.species.names, ...
        complexes_ind, CRN.reactions.details_ind);
    CRN.matrix.S(CRN.species.is_constant==1, :) = 0;  
    
    %% Step 6. Compute conservation laws
    CRN.matrix.Nl = f_compute_semipositive_conservations(CRN.matrix.S);

end