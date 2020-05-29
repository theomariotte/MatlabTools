function freq = freqoct(noct,fmin,fmax,side_freq_fl)
%% 1/n octave band decomposition
% Input param :
%   - noct : decomposition factor (noct = 1 => 1/1 octave bands, noct = 3
%   => 1/3 octave band etc...
%   - fmin : minimum frequency
%   - fmax : maximum frequency
%   - side_freq_fl : include band limit frequencies in freq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fref = 1000.;
dk = 1:15;
lowfact = 2.^(-dk./noct);
upfact = 2.^(dk./noct);

freq_inf = fref .* lowfact;
freq_sup = fref .* upfact;


freq_inf = freq_inf(freq_inf >= fmin);
freq_sup = freq_sup(freq_sup <= fmax);

freq = [freq_inf(end:-1:1) fref freq_sup];


end
