%freqoct_test
clear; clc; close all;

fmin = 60;
fmax = 10000;
n = 1;
fl = 0;

fprintf('No limit bands : \n');
f = freqoct(n,fmin,fmax,fl);
for k = 1 : length(f)
   fprintf('Band %d : f = %6.4e\n',k,f(k)); 
end

% include limit frequencies
fl = 1;
f = freqoct(n,fmin,fmax,fl);
fprintf('Include limit bands : \n');
for k = 1 : length(f)
   fprintf('Band %d | f_inf = %6.4e | f_c = %6.4e | f_sup = %6.4e \n',...
       k,f(k,1),f(k,2),f(k,3)); 
end