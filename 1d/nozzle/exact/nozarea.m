function [a] = nozarea(x)

global L

% straight diverging nozzle
ain = 1.0512;
aout= 1.75;
L   = 10.0;
a = ain + x*(aout-ain)/L;

% another diverging nozzle
%ain = 1.0512;
%aout= 1.75;
%L   = 10.0;
%dc  = 0.8;
%dd  = 4.0;
%db  = (aout - ain)/(tanh(10.0*dc - dd) - tanh(-dd));
%da  = ain - db*tanh(-dd);
%a = da + db*tanh(dc*x - dd);
