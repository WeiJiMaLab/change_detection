function [L P_C_HAT max_P_C_HAT] = run_model_EG(Data,meanvec,lapsevec,priorvec,samples)
Lookup(:,1) = linspace(0,700.92179,70000)';
Lookup(:,2) = besseli(0,Lookup(:,1));
LookupY = Lookup(:,2);
lookup_spacing = 1/(Lookup(2,1)-Lookup(1,1));
lookupEnd = size(Lookup,1);

% precompute log priors
logprior = log(priorvec./(1-priorvec));

% get set sizes
SetSizes = unique(Data(:,5));

L = zeros(length(SetSizes),length(meanvec),length(lapsevec),length(priorvec));
P_C_HAT = zeros([size(L) size(Data,1)]);

% precompute inverse Fisher information function (inverse J)
k = linspace(0,700.92179,600001);
J = k.*(besseli(1,k)./besseli(0,k));

nSamples = samples;

for m = 1:length(SetSizes)
    
    N = SetSizes(m)
    
    % Get set size and delta information
    curr_N_idx = find(Data(:,5)==N);
    nTrials = length(curr_N_idx);
    phi = pi/90*Data(curr_N_idx,56:(55+N));
    theta = pi/90*Data(curr_N_idx,64:(63+N));
    
    C_hat = Data(curr_N_idx,2);
    
    J_vec = meanvec;
    kappa_vec = interp1(J',k',J_vec)';
    kappa_vec_sq = kappa_vec.^2;
    kappa_vec_bessel = besseli(0,kappa_vec);
    
    % loop over mean J
    for jind = 1:length(meanvec)
        
        tic;
        
        % draw noise
        noise1 = reshape(circ_vmrnd(0,kappa_vec(jind),nSamples*N),[],N);
        noise2 = reshape(circ_vmrnd(0,kappa_vec(jind),nSamples*N),[],N);
        
        for i=1:nTrials
            
            % Shift to recycle noise
            noise1 = circshift(noise1,1);
            
            % compute internal representations
            x = bsxfun(@plus,phi(i,:),noise1);
            y = bsxfun(@plus,theta(i,:),noise2);
            
            K3 = sqrt(2*kappa_vec_sq(jind)*(1 + cos(x-y)));
            K3(K3>Lookup(lookupEnd,1)) = Lookup(lookupEnd,1);
            K3(K3<Lookup(1,1)) = Lookup(1,1);
            
            % compute decision variable
            d = log(sum(kappa_vec_bessel(jind)^2./myBessel(K3,lookup_spacing,LookupY),2)) - log(N);                                                                                                %
            
            % add prior
            d = bsxfun(@plus,d,logprior);
            
            % compute p(C_hat = 1 | x,y)
            pC_hat = mean(d>0,1);
            pC_hat = bsxfun(@plus,bsxfun(@times,1-lapsevec',pC_hat),.5*lapsevec');
            
            % compute p(resp | trial info)
            if C_hat(i)==0
                p_resp = 1 - pC_hat;
            else
                p_resp = pC_hat;
            end
            
            p_resp(p_resp==0) = 1/(10*nSamples);
            
            L(m,jind,:,:) = squeeze(L(m,jind,:,:))+log(p_resp);
            
            P_C_HAT(m,jind,:,:,curr_N_idx(i)) = pC_hat;
            
        end
        
        fprintf('N=%d, J1=%2.2f, ETL=%2.2f minutes\n',N,meanvec(jind),(length(meanvec)-jind)*toc/60);

    end
    
end

% print info to screen
max_P_C_HAT = Data;
fprintf('Done!\n',x);