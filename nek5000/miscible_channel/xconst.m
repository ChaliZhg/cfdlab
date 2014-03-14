function [] = xconst(ny1)

system('rm -f output/xconst.dat')
system('cat output/p*.dat >> output/xconst.dat')

data = load('xconst.dat');

% sort based on x coordinate which is in first column
data = sort(data,1,'ascend');

[x,ia,ic] = unique(data(:,1));

% These must be equal
d = ia(2:end) - ia(1:end-1)

m = d(1);
% m/ny1 = no. of elem along x=const

for i=1:length(x)
   jbeg = ia(i);
   if i<length(x)
      jend = ia(i+1)-1; 
      data1 = data(jbeg:jend,:);
   else
      data1 = data(jbeg:end,:);
   end
   data2 = sort(data1,2,'ascend');
end
