function  ris = f_NLPC_restart(x_0, rate_constants, S, Nl, rho, idx_basic_species, ...
    v, ind_one, max_counter, proj)

% Function implementing the Non-Linearly Projected Combined (NLPC)
% algorithm.

% The function 'f_NLPC_restart' takes the following inputs:
% - 'x0' is the starting point for the algorithm
% - 'rate_constants', 'S' 'Nl', 'rho', 'idx_basic_species', 'ind_one', 'v'
% are some data about the stoichiometric surface we're working on and 
% cell's features (in a physiological or mutated state)
% - 'max_counter' is an integer that indicates how many iterations we want
% the algorithm to do every time we choose a new starting point
% - 'proj' is a binry variabile indicating whether to use the classical
% orthogonal projector (proj=1) or the novel non-linear proojector
% (proj=0).

% The function returns the equilibrium computed through the NLPC algorithm
% (with a maximum number of iterations equal to max_counter) having x0 as 
% initial condition and working with the Chemical Rection Network
% defined by rate_constants, S, Nl, rho, idx_basic_species, v, ind_one.

%% Step 1. Define additional parameters within NLPC
toll_cond_init_point = 10^17;
tol = 1e-12;
poss_alpha_2 = ones(1,40);
poss_alpha_2(2) = 0.79;
for i=3:length(poss_alpha_2)
    poss_alpha_2(i) = poss_alpha_2(i-1) * poss_alpha_2(2);
end
poss_alpha = poss_alpha_2(1:20); 
poss_alpha_2 = poss_alpha_2 * 1e3; 
sigma = 10^-4;
sigma_2 = 10^-4;
rho_newcond = 10^-2;
FLAG = 0;

num_try = 500;

ir = 1;

%% Step 2. Compute 'symbolic' form of the Jacobian matrix of v(x; k)
n_species = size(S, 1);
jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);

%% Step 3. Run the algorithm 
while ir < num_try
   
    fprintf('@@@@@@@@@@@  RESTART num %d @@@@@@@@@@ \n', ir)
%%      3.a. Define initial point
%   We draw x_0 on the given stoichiometric compatibility class until
%   cond(J_F(x_0)) is small enough
    if ir > 1
    disp('Initialization')
    aux_cond = Inf; 
    while aux_cond > toll_cond_init_point
        x_0 = f_draw_from_ssurf(Nl, rho, idx_basic_species, [-3, 3]);
        aux_cond = cond(f_evaluate_jacobian(rate_constants, x_0, ...
                    S, idx_basic_species, jacobian_v, Nl));
    end
    disp('Let''s start')
    end
    
%%      3.b. Initialize storing variables
    norm_F_x_nm = zeros(max_counter, 1);
    step_lengths = zeros(max_counter, 1);
    det_F_x = zeros(max_counter, 1);

%%      3.c. External while
    x =  x_0; x(x<0) = 0;
    xnew = x; counter = 0; 
    F_x_new = f_evaluate_mim(rate_constants, xnew, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
    norm_F_x_new = norm(F_x_new);

    while (norm_F_x_new > tol) && (counter < max_counter)
    
        counter = counter + 1;
        x = xnew; F_x = F_x_new; norm_F_x = norm_F_x_new;
        J_x = f_evaluate_jacobian(rate_constants, x, ...
                    S, idx_basic_species, jacobian_v, Nl);

% ************************* Projected Newton ******************************  
        if FLAG == 0 % Newton

            delta = -J_x \ F_x;
            ia = 1;
            while ia <= numel(poss_alpha) 
                alpha = poss_alpha(ia);
                xnew = x + alpha * delta;
                if proj == 0
                    xnew(xnew<0) = x(xnew<0);
                else
                    xnew(xnew<0) = 0;
                end
                F_x_new = f_evaluate_mim(rate_constants, xnew, ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
                norm_F_x_new = norm(F_x_new);
                if norm_F_x_new <= sqrt(1-alpha*sigma)*norm_F_x
                    ia = numel(poss_alpha)+1;
                    flag_save(counter) = 0;
                    FLAG = 0;        
                else
                    FLAG = 1;
                    ia = ia+1;
                end
            end

            if (ia == numel(poss_alpha)+1) && (FLAG == 1)
                xnew = x; F_x_new = F_x; norm_F_x_new = norm_F_x; counter = counter-1;
            else
                % Store some informations
                norm_F_x_nm(counter) = norm_F_x_new;
                step_lengths(counter) = alpha;
                det_F_x(counter) = det(J_x);
%                 fprintf('Iteration %d - f(x) = %2.3e  \n', ...
%                     counter, norm_F_x_nm(counter));
                
                % Condition number of the Jacobian Matrix
                cond_number(counter) = cond(J_x); 
                rcond_number(counter) = rcond(J_x); 
                upd_components(counter) = n_species - sum(xnew==x);
                zeri(counter) = sum(xnew==0);

            end
                
                
% ******************* Projected Gradient Descent **************************
        else 
        
%              disp('********************************************************************************')
            
          
            flag_save(counter) = 1;
            
            delta = - J_x' * F_x;  % delta = - gradiente
            delta_vers = delta / norm(delta);
            ia = 1;

            P_x_grad = x + delta_vers; 
            P_x_grad(P_x_grad < 0) = 0;
            diff_P_x = abs(x - P_x_grad); %
            while ia <= numel(poss_alpha_2)
                alpha = poss_alpha_2(ia);
                xnew = x + alpha * delta_vers;
                diff_P_x_M = diff_P_x(xnew >= 0); %
                diff_P_x_N = diff_P_x(xnew < 0); %
                
                if proj == 0
                    xnew(xnew<0) = x(xnew<0);
                else
                    xnew(xnew<0) = 0;
                end
                F_x_new = f_evaluate_mim(rate_constants, xnew, ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
                norm_F_x_new = norm(F_x_new);
                theta_x = 0.5 * norm_F_x^2;
                theta_x_new = 0.5 * norm_F_x_new^2;
                
                if ((theta_x_new <= theta_x + sigma_2 * (-delta)' * (xnew - x)) ...
                        && (norm(diff_P_x_M) >= rho_newcond * norm(diff_P_x_N)))
                    ia = numel(poss_alpha_2)+1;
                    FLAG = 0;
                else
                    ia = ia+1;
                end
            end
        
        % Store some informations
        norm_F_x_nm(counter) = norm_F_x_new;
        step_lengths(counter) = alpha;
        det_F_x(counter) = det(J_x);
           %fprintf('Iteration %d - f(x) = %2.3e  ia = %d alpha = %2.3e \n', ...
            %   counter, norm_F_x_nm(counter), is, alpha);        
        
        % Condition number of the Jacobian Matrix
        zeri(counter) = sum(xnew==0);
        cond_number(counter) = cond(J_x); 
        rcond_number(counter) = rcond(J_x); 
        upd_components(counter) = n_species - sum(xnew==x);
        
        end 
%         disp("Iteration n. " + counter);
%         disp(sprintf("Norma di F(x): ||F(x)|| = %d ", norm_F_x_new));
    end

%     figure; 
%     semilogy(cond_number,'-o');
%     figure;
%     plot(upd_components);

% 5.c. Store results (for plotting)
x_res = xnew;

% Restart if convergence hasn't been reached
if (counter == max_counter) && (norm_F_x_new > tol)
    ris.cond_number(ir).n = cond_number;
    ris.rcond_number(ir).n = rcond_number;
    ris.upd_components(ir).n = upd_components;
    ris.zeri(ir).n = zeri;
    ris.flag(ir).n = flag_save;
    ir = ir+1; 
    FLAG = 0;
    clear upd_components cond_number zeri flag_save
else
    % Step 6. Store results over run
    ris.x0 = x_0;
    ris.x  = x_res;
    ris.num_iter = counter;
    ris.step_lengths = step_lengths(1:counter);
    ris.norm_F = norm_F_x_nm(counter);
    ris.num_trials = ir;
    ris.upd_components(ir).n = upd_components;
    ris.cond_number(ir).n = cond_number;
    ris.rcond_number(ir).n = rcond_number;
    ris.zeri(ir).n = zeri;
    ris.flag(ir).n = flag_save;
    ir = num_try + 1;
    clear upd_components cond_number zeri rcond_number flag_save
end


end

