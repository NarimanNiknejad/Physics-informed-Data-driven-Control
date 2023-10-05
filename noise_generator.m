function noise = noise_generator(mu,sigma,nStates,nSample,nRealization)
    
    R = chol(sigma);
    
    noise = zeros(nSample,nRealization*nStates);
    
    for ii = 1:nRealization
        
        noise(:,nStates*ii-(nStates-1):nStates*ii) = repmat(mu,nSample,1) + randn(nSample,nStates)*R;
    
    end

end