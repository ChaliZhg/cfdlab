function minmax(U,ubar)
globals;

[N,p] = size(U);
p = p - 1;

for j=1:N
   ff = bezier(U(j,:), xx);
   ucmax(j) = max(ff);
   ucmin(j) = min(ff);
end

%ucmax = max([U(:,1)'; U(:,p+1)'; ubar]);
%ucmin = min([U(:,1)'; U(:,p+1)'; ubar]);

uu = [ucmax(1:N-2); ucmax(2:N-1); ucmax(3:N)];
umax(2:N-1) = max(uu);

uu = [ucmin(1:N-2); ucmin(2:N-1); ucmin(3:N)];
umin(2:N-1) = min(uu);

if periodic==yes
   umin(1) = min([ucmin(N), ucmin(1), ucmin(2)]);
   umin(N) = min([ucmin(N-1), ucmin(N), ucmin(1)]);
   umax(1) = max([ucmax(N), ucmax(1), ucmax(2)]);
   umax(N) = max([ucmax(N-1), ucmax(N), ucmax(1)]);
else
   umin(1) = min([ucmin(1), ucmin(2)]);
   umin(N) = min([ucmin(N-1), ucmin(N)]);
   umax(1) = max([ucmax(1), ucmax(2)]);
   umax(N) = max([ucmax(N-1), ucmax(N)]);
end

%N = length(ubar);
%uu = [ubar(1:N-2); ubar(2:N-1); ubar(3:N)];
%umin(2:N-1) = min(uu);
%umax(2:N-1) = max(uu);
%if periodic==yes
%   umin(1) = min([ubar(N), ubar(1), ubar(2)]);
%   umin(N) = min([ubar(N-1), ubar(N), ubar(1)]);
%   umax(1) = max([ubar(N), ubar(1), ubar(2)]);
%   umax(N) = max([ubar(N-1), ubar(N), ubar(1)]);
%else
%   umin(1) = min([ubar(1), ubar(2)]);
%   umin(N) = min([ubar(N-1), ubar(N)]);
%   umax(1) = max([ubar(1), ubar(2)]);
%   umax(N) = max([ubar(N-1), ubar(N)]);
%end
