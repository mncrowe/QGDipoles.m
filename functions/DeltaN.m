function D = DeltaN(k,l,R,betaU)
% Defines the differential operator DeltaN in Fourier space used for:
%       DeltaN*psi = -(U/a)*sin(theta)*sum[a_j*R_j(r/a)]
%
% (k,l): wavenumbers in (x,y) directions, created using 'create_domain.m'
% R: vector (R_1, R_1, ..., R_N)
% betaU: vector (beta_1/U, beta_2/U, ..., beta_N/U)

if length(R) == 1
    D = -reshape((eps+k.^2+l.^2+1/R^2+betaU),[1 1 length(k) length(l)]);
else
    D = -diag([R(1)^-2 2*reshape(R(2:end-1).^-2,1,[]) R(end)^-2]) - diag(betaU) + ...
        + diag(R(1:end-1).^-2,1) + diag(R(2:end).^-2,-1);
    
    D = repmat(D,1,1,length(k),length(l)) - ...
        repmat(eye(length(R)),1,1,length(k),length(l)).*reshape((eps+k.^2+l.^2),1,1,length(k),length(l));
end

end