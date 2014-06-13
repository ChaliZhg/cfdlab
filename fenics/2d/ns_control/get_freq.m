load history.dat;
t = history(:,1);
ke = history(:,4);
te = history(:,5);

ke = sqrt(ke);
te = sqrt(te);

% skip initial data since there may be some transients
x = ke(200:end);

% make length(x) to be multiple of two
if mod(length(x),2) ~= 0
   x = x(1:end-1);
end

Fs = 1.0/(t(2)-t(1))
x = detrend(x,0);
xdft = fft(x);
freq = 0:Fs/length(x):Fs/2;
xdft = xdft(1:length(x)/2+1);
plot(freq,abs(xdft));
grid on
[~,I] = max(abs(xdft));
fprintf('Maximum occurs at %d Hz.\n',freq(I));
