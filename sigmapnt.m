%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function name: sigmapnt
% Returns the sigma points calculated around the given mean. Matrix
% dimensionality:dX1X2*d+1
% Input: 
% mu=Mean of the state (dX1) vector
% covaria= Covariance of the state dXd matrix
% d = dimension of the state vector

% alpha=Factor deciding the spread of sigma points around the mean
% k= parameter that changes with the distribution of the data.k=2 for
% Gaussian distribution
% lambda = controls the spread of sigma points based on different distributions
% and dimensionality of state vectors
% Output:
% sigmapoints = sigma points around the mean for the particular instance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigmapoints,S] = sigmapnt(mu,covaria,d)
sigmapoints=zeros(d,1,2*d+1);
sigmapoints(:,1,1)=mu;
alpha=0.0001;
k=0;
lambda=alpha^2*(d+k)-d;
gamma = (d + lambda);

%[v,l] = eig((2*d+lambda)*covaria);
S = sqrt(gamma) * sqrtm(covaria);
%shifts = v*sqrt(l)*v^-1;

for i=2:d+1
    sigmapoints(:,1,i)=mu+S(:,(i-1));
end
for i=d+2:2*d+1
    sigmapoints(:,1,i)=mu-S(:,(i-1)-d);
end
    sigmapoints = real(sigmapoints);

if any(~isreal(sigmapoints))
    sigmapoints = real(sigmapoints);
    disp('imaginary sigma points!!!');
end

end

