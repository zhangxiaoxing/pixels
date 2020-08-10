%% AIC grid search is skipped in this implementation 
% fixed history of 2ms *5 window 
%%

ht=10;

% Dimension of input data (L: length, N: number of neurons)
[L,N] = size(X);

% To fit GLM models with different history orders
for suid = 1:N                            % neuron
    % history, when W=2ms
    [bhat{suid}] = glmwin(X,suid,ht,200,2);
end

% To select a model order, calculate AIC
for suid = 1:N
    LLK(suid) = log_likelihood_win(bhat{suid},X,ht,suid,2); % Log-likelihood
end


% Re-optimizing a model after excluding a trigger neuron's effect and then
% Estimating causality matrices based on the likelihood ratio
for target = 1:N
    LLK0(target) = LLK(target);              % Likelihood of full model
    % LLK0(target) = log_likelihood_win(bhat{ht(target),target},X,ht(target),target);
    for trigger = 1:N
        % MLE after excluding trigger neuron
        [bhatc{target,trigger}] = glmcausal(X,target,trigger,ht,200,2);
        
        % Log likelihood obtained using a new GLM parameter and data, which exclude trigger
        LLKC(target,trigger) = log_likelihood_causal(bhatc{target,trigger},X,trigger,ht,target,2);
        
        % Log likelihood ratio
        LLKR(target,trigger) = LLKC(target,trigger) - LLK0(target);
        
        % Sign (excitation and inhibition) of interaction from trigger to target
        % Averaged influence of the spiking history of trigger on target
        SGN(target,trigger) = sign(sum(bhat{target}(ht/2*(trigger-1)+2:ht/2*trigger+1)));
    end 
end

% Granger causality matrix, Phi
Phi = -SGN.*LLKR;

% ==== Significance Test ====
% Causal connectivity matrix, Psi, w/o FDR
D = -2*LLKR;                                     % Deviance difference
alpha = 0.05;
for ichannel = 1:N
    temp1(ichannel,:) = D(ichannel,:) > chi2inv(1-alpha,ht/2);
end
Psi1 = SGN.*temp1;
% keyboard
% Causal connectivity matrix, Psi, w/ FDR
fdrv = 0.05;
temp2 = FDR(D,fdrv,ht);
Psi2 = SGN.*temp2;

% Plot the results
colormap('jet');
figure(1);imagesc(Phi,[-std(Phi(:)),std(Phi(:))]);xlabel('Triggers');ylabel('Targets');
colormap('jet');
figure(2);imagesc(Psi1,[-1,1]);xlabel('Triggers');ylabel('Targets');
colormap('jet');
figure(3);imagesc(Psi2,[-1,1]);xlabel('Triggers');ylabel('Targets');

% Save results
save('CausalMaps','bhatc','LLK0','LLKC','LLKR','D','SGN','Phi','Psi1','Psi2');