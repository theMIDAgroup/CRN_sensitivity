function [SSI_k, SSI_c] = f_compute_SSI(idx_sp, x_eq, rate_constants, S, Nl, rho, idx_basic_species, v)

%Function computing the statistical sensitivity indices with respect to
% - reactions
% - conservation laws

% The function 'f_compute_SSI' takes the following inputs:
% - 'idx_sp' is the set of selected species to compute the SSI
% - 'x_eq' is the equilibrium point of the network
% - 'rate_constants', 'S' 'Nl', 'rho', 'idx_basic_species', 'v'
% are some data about the stoichiometric surface we're working on and 
% cell's features (in a physiological or mutated state)
% - 'mut_lof' is a boolean determining if the network underwent a loss of
% function mutation: 1 yes, 0 no

% The function returns the SSI computed through the having x_eq as 
% the equilibrium point for a specific set of species (determined by idx_sp)
% and working with the Chemical Rection Network defined by rate_constants, 
% S, Nl, rho, idx_basic_species, and v.

% Specifically
% SSI_k : statistical sensitivity indices with respect to reactions
% SSI_c : statistical sensitivity indices with respect to conservation laws


% 1. check x_eq has non-zero components
if isempty(find(x_eq(idx_sp)==0, 1))==0
    error('The equilibrium point has 0 components')
else
    % 2. compute relative sensitivity matrix
    % 2.1 general conditions

    n_species=numel(x_eq);
    ind_one=n_species+1;
    n_cons_laws=size(Nl,1);

    % 2.2 computation of general sensitivity matrix

    jacobian_v= f_compute_analytic_jacobian_v(v, n_species, ind_one);
    
    eval_jac = f_evaluate_jacobian_neworder(rate_constants, x_eq, ...
        S, idx_basic_species, jacobian_v, Nl);
    
    inv_eval_jac=inv(eval_jac);

    % 2.3. check inversion of jacobian

    if inv_eval_jac==Inf
        error('Cannot continue computing... ending function')
    else
        
        % 2.4. computation of sensitivity matrix

        eval_jac_k=f_evaluate_jacobian_k(v, x_eq, S, idx_basic_species);
        eval_jac_c=f_evaluate_jacobian_c(n_species, n_cons_laws);
        
        size(eval_jac_c)
        
        Jac_k_xeq= -inv_eval_jac*eval_jac_k;
        Jac_c_xeq= -inv_eval_jac*eval_jac_c;
        
        % 2.5. row (=species ) selection
        Jac_k_xeq=Jac_k_xeq(idx_sp, :);
        Jac_c_xeq=Jac_c_xeq(idx_sp, :);
        x_eq=x_eq(idx_sp);
       
        % 2.6. computation of relative sensitivity matrix
        n_species_new=numel(x_eq);
        selection_k=min(n_species_new, numel(rate_constants));
        selection_c=min(n_species_new, n_cons_laws);
        
        L_k=(1./x_eq).*(Jac_k_xeq).*(rate_constants');
        [~,D_k,V_k]=svd(L_k);
   
        eing_k=diag(D_k(1:selection_k, 1:selection_k)).^2;
        sum_eing_k=sum(eing_k);
        SSI_k =  ((V_k(:,1:selection_k).^2)*(eing_k))/sum_eing_k; %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
        
        size(Jac_c_xeq)
        
        L_c=(1./x_eq).*(Jac_c_xeq).*(rho');
        [~,D_c,V_c]=svd(L_c);
        
        eing_c=diag(D_c(1:selection_c, 1:selection_c)).^2;
        sum_eing_c=sum(eing_c);
        
        SSI_c =  ((V_c(:,1:selection_c).^2)*(eing_c))/sum_eing_c; %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
    
    end
end


