close all
clear all

nx=5;
ny=5;

figure(1)
for i=1:nx
   for j=1:ny
      plot([i],[j],'o')
      hold on
   end
end

for i=0:nx
   j=0
   plot([i+1/2], [j+1/2], '*')
end
