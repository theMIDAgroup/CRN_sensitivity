function  ris = f_PNG_restart(x_0, rate_constants, S, Nl, rho, idx_basic_species, ...
    v, ind_one, max_counter)

%% TODO:
% 1. Il calcolo/salvataggio di alcune variabili potrebbe essere reso 
%    opzionale

% Input:
% - S
% - idx_basic_species
% - x_0

% The function 'newtonGD_restart' takes the following inputs:
% - 'MIM' is a struct that contains all cell's features (in a physiological or mutated state)
% - 'ris_dynamic' is a struct that contains the results obtained working on the
% same cell, but using the dynamic approach
% - 'num_run' is an integer that indicates how many experiments we want to
% do on the same stoichiometric matrix
% - 'max_counter' is an integer that indicates how many iterations we want
% the algorithm to do every time we choose a new starting point

%% Step 1. Define additional parameters within PNG
toll_cond_init_point = 10^17; %22
tol = 1e-12; 
poss_alpha = logspace(0, -2, 20); %20
poss_alpha_2 = 1e3 * logspace(0, -2, 20);
poss_alpha_3 = 1e1 * logspace(0, -2, 20);
poss_alpha_2 = [poss_alpha_2  poss_alpha_3];
sigma = 10^-4;
sigma_2 = 10^-4;
FLAG = 0;

num_try = 100;

ir = 1;

%% Step 2. Compute 'symbolic' form of the Jacobian matrix of v(x; k)
n_species = size(S, 1);
jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);

%% Step 4. Run the algorithm 
while ir < num_try
   
    fprintf('@@@@@@@@@@@  RUN num %d @@@@@@@@@@ \n', ir)
%%      4.a. Define initial point
%   We draw x_0 on the given stoichiometric compatibility class until
%   cond(J_F(x_0)) is small enough
    if ir >= 1
    disp('Initialization')
    aux_cond = Inf; 
    while aux_cond > toll_cond_init_point
        x_0 = f_draw_from_ssurf(Nl, rho, idx_basic_species, [-3, 3]);
        aux_cond = cond(f_evaluate_jacobian(rate_constants, x_0, ...
                    S, idx_basic_species, jacobian_v, Nl));
    end
    disp('Let''s start')
    end
    
%%      4.b. Initialize storing variables
    norm_F_x_nm = zeros(max_counter, 1);
    step_lengths = zeros(max_counter, 1);
    det_F_x = zeros(max_counter, 1);
    norm_grad = zeros(1,max_counter);

%%      4.c. Run newton method 
    x =  x_0; x(x<0) = 0;
    xnew = x; counter = 0; 
    F_x_new = f_evaluate_mim(rate_constants, xnew, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
    norm_F_x_new = norm(F_x_new);
    
    j = 1;

    while (norm_F_x_new > tol) && (counter <= max_counter)
    
        counter = counter + 1;
        x = xnew; F_x = F_x_new; norm_F_x = norm_F_x_new;
%new_order?
        J_x = f_evaluate_jacobian(rate_constants, x, ...
                    S, idx_basic_species, jacobian_v, Nl);

% ************************* Projected Newton ******************************  
        if FLAG == 0 % Newton

            delta = -J_x \ F_x;
            ia = 1;
            while ia <= numel(poss_alpha) 
                alpha = poss_alpha(ia);
                xnew = x + alpha * delta;
                xnew(xnew<0) = x(xnew<0);
                F_x_new = f_evaluate_mim(rate_constants, xnew, ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
                norm_F_x_new = norm(F_x_new);
                if norm_F_x_new <= sqrt(1-alpha*sigma)*norm_F_x
                    ia = numel(poss_alpha)+1;
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
            end

% ******************* Projected Gradient Descent **************************
        else 
        
%             disp('********************************************************************************')

            delta = - J_x' * F_x;
            delta_vers = delta / norm(delta);
            ia = 1;
            while ia < numel(poss_alpha_2)
                alpha = poss_alpha_2(ia);
                xnew = x + alpha * delta_vers;
                xnew(xnew<0) = x(xnew<0);
                F_x_new = f_evaluate_mim(rate_constants, xnew, ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
                norm_F_x_new = norm(F_x_new);
                theta_x = 0.5 * norm_F_x^2;
                theta_x_new = 0.5 * norm_F_x_new^2;
                if theta_x_new <= theta_x + sigma_2 * (-delta_vers)' * (xnew - x)
                    is = ia;
                    ia = numel(poss_alpha_2);
                    FLAG = 0;
                    norm_grad(j) = norm(delta);
                    j = j+1;
                else
                    ia = ia+1;
                    alpha = poss_alpha_2(ia);
                    is = ia;
                end
            end

        % Store some informations
        norm_F_x_nm(counter) = norm_F_x_new;
        step_lengths(counter) = alpha;
        det_F_x(counter) = det(J_x);
%         fprintf('Iteration %d - f(x) = %2.3e  ia = %d \n', ...
%             counter, norm_F_x_nm(counter), is);        
        
        end 
    end

% 5.c. Store results (for plotting)
x_res = xnew;

% Restart if convergence hasn't been reached
if (counter == max_counter+1) && (norm_F_x_new > tol)
    ir = ir+1; 
else
    % Step 6. Store results over run
    ris.x0 = x_0;
    ris.x  = x_res;
    ris.num_iter = counter;
    ris.step_lengths = step_lengths(1:counter);
    ris.norm_F = norm_F_x_nm(counter);
    ris.num_trials = ir;
    ir = num_try + 1;
end


end

% struct.ris = ris;
% 
% %% Check: how many times dynamic and Newton-GD get to the same result?
% 
% 
% %% Norm of F in the equilibrium point
% 
% mean_norm_F = 0;
% 
% for i = 1:num_run
%     mean_norm_F = ris(i).norm_F + mean_norm_F;
% end
% 
% struct.all.mean_norm_F = mean_norm_F / num_run;
% 
% %% Convergence rate
% 
% num_trials = 0;
% 
% for i=1:num_run
%     num_trials = num_trials + ris(i).num_trials;
% end
% 
% struct.all.effect_conv = num_run/num_trials;
