function freq = freqoct(noct,fmin,fmax,side_freq_fl)
%% 1/n octave band decomposition
% Input param :
%   - noct : decomposition factor (noct = 1 => 1/1 octave bands, noct = 3
%   => 1/3 octave band etc...
%   - fmin : minimum frequency
%   - fmax : maximum frequency
%   - side_freq_fl : include band limit frequencies in freq (usefull for
%   band filter design)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4, side_freq_fl = []; end
if nargin < 3, fmax = []; end
if nargin < 2, fmin = []; end
if nargin < 1, noct = []; end

% default : 1/1 octave band decomposition between 63 Hz and 8000 Hz without
% limit frequencies
if isempty(side_freq_fl), side_freq_fl = 0; end
if isempty(fmax), fmax = 8000.; end
if isempty(fmin), fmin = 62.; end
if isempty(noct), noct = 1; end

fref = 1000.;
% TODO : this can be upgraded to be automated
dk = 1:16;
lowfact = 2.^(-dk./noct);
upfact = 2.^(dk./noct);

freq_inf = fref .* lowfact;
freq_sup = fref .* upfact;

freq_inf = freq_inf(freq_inf >= fmin);
freq_sup = freq_sup(freq_sup <= fmax);

freq = [freq_inf(end:-1:1) fref freq_sup];
freq = shiftdim(freq);

if side_freq_fl
   
    freq_inf = freq.* 2^(-1/2/noct);
    freq_sup = freq.* 2^(+1/2/noct);
    freq = [freq_inf freq freq_sup];
    
end


end
