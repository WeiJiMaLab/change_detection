function [L P_C_HAT max_P_C_HAT] = run_model_SR(Data,meanvec,kvec,priorvec,samples)

Lookup = linspace(0,5000,5000000)';
LookupY = besseli(0,Lookup(:,1),1);
lookup_spacing = 1/(Lookup(2)- Lookup(1));

% get priors
logprior = log(priorvec./(1-priorvec));

L = zeros(length(meanvec),length(kvec),length(priorvec));
P_C_HAT = zeros([size(L) size(Data,1)]);

% precompute inverse Fisher information function (inverse J)
k = linspace(0,5000,10000000)';
J = k.*(besseli(1,k,1)./besseli(0,k,1));

nSamples = samples;

% get set sizes
SetSizes = unique(Data(:,5));

all_kappa_mat = interp1(J,k,bsxfun(@rdivide,meanvec',kvec));
all_kappa_bessel = besseli(0,all_kappa_mat,1);

for m = 1:length(SetSizes)
    
    N = SetSizes(m);
    
    % Initialize priors and Likelihood matrix
    curr_N_idx = find(Data(:,5)==N);
    nTrials = length(curr_N_idx);
    phi = pi/90*Data(curr_N_idx,56:(55+N));
    theta = pi/90*Data(curr_N_idx,64:(63+N));
    nancount = zeros(1,length(meanvec));

    % loop over mean J
    for jind = 1:length(meanvec)
        
        tic;
        
        % loop over theta parameter
        for kind = 1:length(kvec)
            
            K = kvec(kind);
            
            if K>N
                K=N;
            end
            
            kappa = all_kappa_mat(jind,K);
            kappa_sq = kappa^2;
            kappa_bessel = all_kappa_bessel(jind,K);
            kappa_bessel_log = log(kappa_bessel);
            
            
            noise1 = reshape(circ_vmrnd(0,kappa,nSamples*N),nSamples,N);
            noise2 = reshape(circ_vmrnd(0,kappa,nSamples*N),nSamples,N);
            
            [Y I] = sort(rand(nSamples,N),2,'ascend');
            keep_mask = I <= K;
            
            for i=1:nTrials
                
                noise1 = noise1([2:end 1],:);
                keep_mask = keep_mask([2:end 1],:);
                
                costerm = cos(bsxfun(@plus,phi(i,:)-theta(i,:),noise1-noise2));
                K3 = sqrt(2*kappa_sq*(1+costerm));
                if ~isreal(K3)
                    fprintf('\nOMGimaginary!\n');
                end
                
                K3_bessel = myBessel(K3,lookup_spacing,LookupY);
                d = log(row_sum(exp(2*kappa_bessel_log-log(K3_bessel)+2*kappa-K3).*keep_mask))-log(K);
                
                d = bsxfun(@plus,d,logprior);
                
                % compute p(C_hat = 1 | x,y)
                pC_hat = sum(d>0)/nSamples;
                
                P_C_HAT(jind,kind,:,curr_N_idx(i)) = shiftdim(pC_hat,-1);
                
            end
        end
        
        fprintf('N=%d, J1=%2.2f, ETL=%2.2f minutes\n',N,meanvec(jind),(length(meanvec)-jind)*toc/60);
          
    end
    
end

% print info to screen
max_P_C_HAT = Data; % disregard
fprintf('Done!\n');