% Hautus Criterion

function [nu] = hautus(A,B,M)

n = length(A(:,1));

i = 1;
[V,D] = eigs(A',M',i,'la');
lambda = max(diag(D));

% In case the system is already stable
if (lambda < 0)
    disp('System is trivially stabilizable as all eigenvalues are stable')
    nu = 0;
    return 
end
    
% Determining the unstable eigenvalues and eigenvectors    
while (lambda >= 0 && i <=n)
   i = i+1;
   [V,D] = eigs(A',M',i,'la');
   lambda = min(diag(D));
end

% Incase the while loop exited because lambda < 0
if (lambda < 0)
   V = V(:,1:i-1);
   D = D(1:i-1,1:i-1);
end
    
tol = 1.0e-12;

count = 0;
nu = length(V(1,:)); %number of unstable eigenvectors

for j = 1:nu
   if( norm(B'*V(:,j)) > tol)
       count = count + 1;
   end
end

if count == nu
   disp('system is stabilizable')
else
   disp('system is not stabilizable')
end

