function [L P_C_HAT max_P_C_HAT] = run_model_IL(Data,lapsevec,kvec,qvec)

L = zeros(length(kvec),length(lapsevec),length(qvec));
P_C_HAT = zeros([size(L) size(Data,1)]);

% get set sizes
SetSizes = unique(Data(:,5));
ss_idx = 0;

for m = 1:length(SetSizes)
    
    N = SetSizes(m);
    
    % Initialize priors and Likelihood matrix
    curr_N_idx = find(Data(:,5)==N);
    nTrials = length(curr_N_idx);
    phi = pi/90*Data(curr_N_idx,56:(55+N));
    theta = pi/90*Data(curr_N_idx,64:(63+N));

    Delta = abs(sum(circ_dist(phi,theta),2)) > 1e-10;
    
    C_hat = Data(curr_N_idx,2);
    
    for qind = 1:length(qvec)
        
        q = qvec(qind);
        
        
        % loop over K parameter
        for kind = 1:length(kvec)
            
            K = kvec(kind);
            
            % loop over lapse rate parameter
            for lind = 1:length(lapsevec)
                
                lapse = lapsevec(lind);
                for i=1:nTrials
                    
                    % compute p(resp | trial info)
                    if K>N
                        K=N;
                    end
                    
                    % Probability of responding "yes"
                    if Delta(i)==1
                        pC_hat = (1-lapse)*K/N + (1-K/N)*q;
                    else
                        pC_hat = q;
                    end

                    % Probability of the subject's response
                    if C_hat(i)==1
                        p_resp = pC_hat;                            
                    else
                        p_resp = 1-pC_hat;
                    end
                    
                    % save p_resp
                    P_C_HAT(kind,lind,qind,curr_N_idx(i)) = pC_hat;
                    
                    L(kind,lind,qind) = L(kind,lind,qind) + log(p_resp);
                end
            end
        end
    end
    ss_idx = ss_idx+nTrials;
end

% print info to screen
max_P_C_HAT = Data; % disregard
fprintf('Done!\n');