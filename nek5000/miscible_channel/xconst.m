% ny1 = number of GLL points in y direction
% ii  = index of x = const line
function [] = xconst(ny1,ii)

system('rm -f output/xconst.dat');
system('cat output/p*.dat >> output/xconst.dat');
pause(2)

data = load('output/xconst.dat');

% sort based on x coordinate which is in first column
data = sortrows(data,1);

%figure(1)
%plot(data(:,1),data(:,2),'-')

[x,ia,ic] = unique(data(:,1),'stable');

% These must be equal
d = ia(2:end) - ia(1:end-1);
for i=1:length(d)-1
   assert(d(i+1)==d(i));
end

m = d(1);
assert(mod(m,ny1)==0);
nelem = m/ny1;
% m/ny1 = no. of elem along x=const

for i=1:length(x)
   %[i,x(i)]
   jbeg = ia(i);
   if i<length(x)
      jend = ia(i+1)-1; 
      data1 = data(jbeg:jend,:);
   else
      data1 = data(jbeg:end,:);
   end
   data2 = sortrows(data1,2);
   % check element ends are correct
   for e=1:nelem-1
      assert(data2(e*ny1,2)==data2(e*ny1+1,2),'%f ~= %f\n', ...
             data2(e*ny1,2), data2(e*ny1+1,2));
      assert(data2(e*ny1,3)==data2(e*ny1+1,3),'%f ~= %f\n', ...
             data2(e*ny1,3), data2(e*ny1+1,3));
   end

%  figure(2)
%  plot(data2(:,3), data2(:,2), 'o-')
%  figure(3)
%  plot(data2(:,6), data2(:,2), '-')
%  pause

   if i==ii
      fprintf(1,'Saving at x = %f into xconst.mat\n', x(i))
      save('xconst','data2')
      break
   end
end

fprintf(1,'i index %d not found\n', ii)
fprintf(1,'Maximum i index is %d\n', length(x))
