% Hautus Criterion
function [nu] = hautus(A,B)

n = length(A(:,1));

i = 1;
[V,D] = eig(A');

% Determining the unstable eigenvalues and eigenvectors    

count = 0;
nu = 0;
nz = 0;
tol = 1e-14;
for j = 1:n
   if(D(j,j) >= 0)
       if(D(i,j) > 0)
           nu = nu + 1;
       else
           nz = nz + 1;
       end
       if( norm(B'*V(:,j)) > tol)
          count = count + 1;
       end
   end   
end

fprintf('Number of positive eigenvalues = %d \n',nu)
fprintf('Number of zero eigenvalues = %d\n\n',nz)

if count == nu + nz
   disp('system is stabilizable')
else
   disp('system is not stabilizable')
end

