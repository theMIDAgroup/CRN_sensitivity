function SSI = f_compute_SSI_tot(idx_sp, x_eq, rate_constants, S, Nl, rho, idx_basic_species, v)

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
        error('Cannot continue computing... stop function')
    else
        
        % 2.4. computation of sensitivity matrix

        eval_jac_k=f_evaluate_jacobian_k(v, x_eq, S, idx_basic_species);
        eval_jac_c=f_evaluate_jacobian_c(n_species, n_cons_laws);
        
        eval_jac_h = [eval_jac_k, eval_jac_c];
        
        Jac_h_xeq= -inv_eval_jac*eval_jac_h;
        
        % 2.5. row (=species ) selection
        Jac_h_xeq=Jac_h_xeq(idx_sp, :);
        
        x_eq=x_eq(idx_sp);
       
        % 2.6. computation of relative sensitivity matrix

        selection=min(numel(idx_sp), (numel(rate_constants)+numel(rho)));
        
        L=(1./x_eq).*(Jac_h_xeq).*([rate_constants; rho]');
        [~,D,V]=svd(L);
        
        lambda=diag(D(1:selection, 1:selection)).^2;
        sum_lambda=sum(lambda);

        SSI =  ((V(:,1:selection).^2)*lambda)/sum_lambda; %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
        
         
    end
end


