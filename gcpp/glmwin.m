function [beta_new] = glmwin(X,n,ht,windowWidth,w)

%================================================================
%                GLM fitting based on submatrices
%================================================================
%
%  This code is made for the case when input matrix X is too large
%  X is partioned into small submatrices of (k x 1)-dimension
%  This code is based on bnlrCG.m (Demba) and 
%
%   References:
%      Dobson, A.J. (1990), "An Introduction to Generalized Linear
%         Models," CRC Press.
%      McCullagh, P., and J.A. Nelder (1990), "Generalized Linear
%         Models," CRC Press.
%
% Input arguments:
%           X: measurement data (# samples x # neurons)
%           n: index number of input (target neuron) to analyze
%          ht: model order (using AIC or BIC)
%           k: one of divisors for K
%           w: duration of non-overlapping spike counting window
%
% Output arguments:
%     beta_new: estimated GLM parameters
%
%================================================================
% Sanggyun Kim
% Neuroscience Statistics Research Lab (BCS MIT)
% April 13. 2009
%================================================================

% Spike counting window
WIN = zeros(ht/w,ht);
for iwin = 1:ht/w
    WIN(iwin,(iwin-1)*w+1:iwin*w) = 1;
end

% CG parameters
cgeps = 1e-3;
cgmax = 30;

% LR parameters
Irmax = 100;
Ireps = 0.05;

% Design matrix, including DC column of all ones (1st or last)
Xnew = X(ht+1:end,:);
[XnewFullWidth,Q] = size(Xnew);% samples, SUs
NxStepPlusBias = Q*ht/w +1;
for currWindow = 1:XnewFullWidth/windowWidth
    for currIdxInWindow = 1:windowWidth - ht
        temp = [1];
        for SUid = 1:Q
            temp = [temp (WIN*X(ht+(currWindow-1)*windowWidth+currIdxInWindow-1:...
            -1:(currWindow-1)*windowWidth+currIdxInWindow,SUid))']; %jj is neuron, transverse length of ht for each neuron
        end
        Xsub{currWindow,1}(currIdxInWindow,1:NxStepPlusBias) = temp;
    end
end

% Making output matrix Ysub{}
for currWindow = 1:XnewFullWidth/windowWidth
    Ysub{currWindow} = Xnew(windowWidth*(currWindow-1)+1:windowWidth*currWindow-ht,n);
end

% Logistic regression
i = 0;
% Initialization
% P = length(Xsub{1,1});
beta_old = zeros(NxStepPlusBias,1);%neuron * time point + inception
for currWindow = 1:XnewFullWidth/windowWidth  %number of small windows
%     eta{currWindow} = zeros(windowWidth-ht,1);
    eta{currWindow} = Xsub{currWindow}*beta_old(1:NxStepPlusBias);
    %% equiv?
    %% eta{kk}=Xsub{kk}*beta_old;
    musub{currWindow} = exp(eta{currWindow})./(1+exp(eta{currWindow}));
    Wsub{currWindow} = diag(musub{currWindow}).*diag(1-musub{currWindow});
    zsub{currWindow} = eta{currWindow} + (Ysub{currWindow}-musub{currWindow}).*(1./diag(Wsub{currWindow}));
end

% Scaled deviance
devold = 0;
for currWindow = 1:XnewFullWidth/windowWidth
    devold = devold - 2*(Ysub{currWindow}'*log(musub{currWindow})+(1-Ysub{currWindow})'*log(1-musub{currWindow}));
end
devnew = 0;
devdiff = abs(devnew - devold);

% Do CG -> beta_new, i.e. solve for beta_new: X'WX*beta_new =
% X'Wz(beta_old) using CG
while (i < Irmax && devdiff > Ireps)
    
    A(1:NxStepPlusBias,1:NxStepPlusBias) = zeros(NxStepPlusBias,NxStepPlusBias);
    for currWindow = 1:XnewFullWidth/windowWidth
        A(1:NxStepPlusBias,1:NxStepPlusBias) = A(1:NxStepPlusBias,1:NxStepPlusBias) + Xsub{currWindow}'*Wsub{currWindow}*Xsub{currWindow};
    end
    
    %A = A + A' - diag(diag(A));

    b(1:NxStepPlusBias,1) = zeros(NxStepPlusBias,1);
    for currWindow = 1:XnewFullWidth/windowWidth
        b(1:NxStepPlusBias,1) = b(1:NxStepPlusBias,1) + Xsub{currWindow}'*Wsub{currWindow}*zsub{currWindow};
    end
    

    % Conjugate gradient method for symmetric postive definite matrix A
    beta_new = cgs(A,b,cgeps,cgmax,[],[],beta_old);
    beta_old = beta_new;

    for currWindow = 1:XnewFullWidth/windowWidth
%         eta{currWindow} = zeros(windowWidth-ht,1);
        eta{currWindow} = Xsub{currWindow}*beta_old(1:NxStepPlusBias);
        musub{currWindow} = exp(eta{currWindow})./(1+exp(eta{currWindow}));
        Wsub{currWindow} = diag(musub{currWindow}).*diag(1-musub{currWindow});
        zsub{currWindow} = eta{currWindow} + (Ysub{currWindow}-musub{currWindow}).*(1./diag(Wsub{currWindow}));
    end

    % Scaled deviance
    devnew = 0;
    for currWindow = 1:XnewFullWidth/windowWidth
        devnew = devnew - 2*(Ysub{currWindow}'*log(musub{currWindow})+(1-Ysub{currWindow})'*log(1-musub{currWindow}));
    end
    devdiff = abs(devnew - devold);
    devold = devnew;
    
    i = i+1;
    
end

% % Compute additional statistics
% stats.dfe = 0;
% stats.s = 0;
% stats.sfit = 0;
% stats.covb = inv(A);
% stats.se = sqrt(diag(stats.covb));
% stats.coeffcorr = stats.covb./sqrt((repmat(diag(stats.covb),1,p).*repmat(diag(stats.covb)',p,1)));
% stats.t = 0;
% stats.p = 0;
% stats.resid = 0;
% stats.residp = 0;
% stats.residd = 0;
% stats.resida = 0;