function [L P_C_HAT max_P_C_HAT] = run_model_SA(Data,meanvec,kvec,priorvec,samples)

% for myBessel
Lookup(:,1) = linspace(0,700.92179,35000)';
Lookup(:,2) = besseli(0,Lookup(:,1));
LookupEnd = Lookup(size(Lookup,1),1);
lookup_spacing = (1/(Lookup(2,1)-Lookup(1,1)));
LookupY = Lookup(:,2);

% get priors
logprior = log(priorvec./(1-priorvec));

% init
L = zeros(length(meanvec),length(kvec),length(priorvec));
P_C_HAT = zeros([size(L) size(Data,1)]);

% precompute inverse Fisher information function (inverse J)
k = linspace(eps,700.92179,6001);
J = k.*(besseli(1,k)./besseli(0,k));

nSamples = samples;

% get set sizes
SetSizes = unique(Data(:,5));

for m = 1:length(SetSizes)
    
    % take the current set size of data
    N = SetSizes(m);
    curr_N_idx = find(Data(:,5)==N);
    
    % get the number of trials and orientations
    nTrials = length(curr_N_idx);
    phi = pi/90*Data(curr_N_idx,56:(55+N));
    theta = pi/90*Data(curr_N_idx,64:(63+N));
    
    % get the responses
    C_hat = Data(curr_N_idx,2);
    
    % loop over mean J
    for jind = 1:length(meanvec)
        
        tic;
        
        J1 = meanvec(jind); % Fisher info for 1 slot
        kappa = interp1(J',k',J1); % kappa for 1 slot
        kappa_bessel = besseli(0,kappa);
        
        % loop over K parameter (number of slots)
        for kind = 1:length(kvec)
            
            % get kappa high and low resource for K>N
            K = kvec(kind);
            kappa_low  = interp1(J',k',floor(K/N)*J1);
            kappa_high = interp1(J',k',(floor(K/N)+1)*J1);
            
            % get besseli for faster processing
            kappa_low_bessel = besseli(0,kappa_low);
            kappa_high_bessel = besseli(0,kappa_high);
            
            % CASE 1: Fewer slots than items
            if K<=N
                
                % create idx matrix which will tell which stimuli to keep
                [dummy,keep_idx] = sort(rand([nSamples N nTrials]),2);
                
                % decide which ones to keep
                idx_temp = keep_idx<= K;
                
                % Create noise matrices for first and second screen
                noise1 = reshape(circ_vmrnd(0,kappa,nSamples*N),nSamples,[]);
                noise2 = reshape(circ_vmrnd(0,kappa,nSamples*N),nSamples,[]);
                
                for i=1:nTrials
                    noise1 = circshift(noise1,1);
                    
                    % compute difference between x and y
                    real_delta = phi(i,:)-theta(i,:);
                    noise_delta = noise1 - noise2;
                    full_delta = bsxfun(@plus,real_delta,noise_delta);
                   
                    % calculate Kappa_c
                    K3 = kappa * sqrt(2 + 2*cos(full_delta));
                    K3(K3>LookupEnd) = LookupEnd;
                    K3(K3<Lookup(1,1)) = Lookup(1,1);
                    
                    % compute decision variable
                    d = log(sum(idx_temp(:,:,i).*(kappa_bessel)^2./myBessel(K3,lookup_spacing,LookupY),2)) - log(K);
                    d = bsxfun(@plus,d,logprior);
                    
                    % compute p(C_hat = 1 | x,y)
                    pC_hat = mean(d>0,1);
                    
                    % compute p(resp | trial info)
                    if C_hat(i)<0.5
                        p_resp = 1 - pC_hat;
                    else
                        p_resp = pC_hat;
                    end
                    
                    % take care of 0-responses
                    p_resp(p_resp==0) = 1/(10*nSamples);
                    
                    % store results
                    L(jind,kind,:) = L(jind,kind,:) + shiftdim(log(p_resp),-1);
                    
                    P_C_HAT(jind,kind,:,curr_N_idx(i)) = shiftdim(pC_hat',-2);
                    
                end
                % CASE 2: At least as many slots as items
            else
                % number of high Kappa stimuli
                nHigh = mod(K,N);
                
                % create idx matrix which will tell which stimuli to keep
                [dummy,keep_idx] = sort(rand([nSamples N nTrials]),2);
                
                % decide which ones to keep
                idx_temp = keep_idx <= nHigh;
                
                % Create noise matrices for first and second screen
                noise1_low = reshape(circ_vmrnd(0,kappa_low,nSamples*N),nSamples,[]);
                noise2_low = reshape(circ_vmrnd(0,kappa_low,nSamples*N),nSamples,[]);
                noise1_high = reshape(circ_vmrnd(0,kappa_high,nSamples*N),nSamples,[]);
                noise2_high = reshape(circ_vmrnd(0,kappa_high,nSamples*N),nSamples,[]);
                
                for i = 1:nTrials
                    noise1_low = circshift(noise1_low,1);
                    noise1_high = circshift(noise1_high,-1);
                    
                    % randomly set low/high (nHigh 1's per row, N-nHigh 0's per row)
                    noise1 = idx_temp(:,:,i).*noise1_high + (1-idx_temp(:,:,i)).*noise1_low;
                    noise2 = idx_temp(:,:,i).*noise2_high + (1-idx_temp(:,:,i)).*noise2_low;
                    
                    % set kappas and their bessel values
                    kappa_mat = kappa_high*idx_temp(:,:,i)+kappa_low*(1-idx_temp(:,:,i));
                    
                    bessel_mat = kappa_high_bessel*idx_temp(:,:,i)+kappa_low_bessel*(1-idx_temp(:,:,i));
                    
                    % compute difference between x and y
                    real_delta = (phi(i,:)-theta(i,:));
                    noise_delta = (noise1-noise2);
                    full_delta = bsxfun(@plus,real_delta,noise_delta);
                    
                    % get decision variable
                    K3 = kappa_mat .* sqrt(2 + 2*cos(full_delta));
                    K3(K3>LookupEnd) = LookupEnd;
                    K3(K3<Lookup(1,1)) = Lookup(1,1);
                    
                    % compute decision variable
                    d = log(sum((bessel_mat.^2)./myBessel(K3,lookup_spacing,LookupY),2)) - log(N);
                    d = bsxfun(@plus,d,logprior);
                    
                    % compute p(C_hat = 1 | x,y)
                    pC_hat = mean(d>0);
                    
                    % compute p(resp | trial info)
                    if C_hat(i)<.5
                        p_resp = 1 - pC_hat;
                    else
                        p_resp = pC_hat;
                    end
                    
                    % take care of 0-responses
                    p_resp(p_resp==0) = 1/(10*nSamples);
                    
                    % store results
                    
                    L(jind,kind,:) = L(jind,kind,:) + shiftdim(log(p_resp),-1);
                    
                    P_C_HAT(jind,kind,:,curr_N_idx(i)) = shiftdim(pC_hat,-1);
                    
                end
                
            end
        end
        
        fprintf('N=%d, J1=%2.2f, ETL=%2.2f minutes\n',N,meanvec(jind),(length(meanvec)-jind)*toc/60);

        
    end
end

% print info to screen
max_P_C_HAT = Data; % disregard
fprintf('Done!\n');