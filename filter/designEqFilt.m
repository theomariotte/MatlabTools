function filterBank = designEqFilt(filtParam)
%% Banc de filtres passe-bande - décomposition en 1/n d'octave
% 
% filtParam : 
%   - 'fmin' : lowest band center frequency
%   - 'fmax' : 'highest band center frequency
%   - 'fs'   : sampling rate
%   - 'noct' : band decomposition (e.g. noct = 3 means 1/3 octave bands)
%   - 'energyfl' : double filtering to conserve energy (useful for
%   equalization)
%   - 'order' : filter order
%   - plotfl : (0) no plot ; (1) plot frequency response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(filtParam,'fmin'), filtParam.fmin = []; end
if ~isfield(filtParam,'fmax'), filtParam.fmax = []; end
if ~isfield(filtParam,'fs'), filtParam.fmin = []; end
if ~isfield(filtParam,'noct'), filtParam.noct = []; end
if ~isfield(filtParam,'order'), filtParam.order = 6; end
if ~isfield(filtParam,'energyfl'), filtParam.energyfl = 0; end
if ~isfield(filtParam,'plotfl'), filtParam.plotfl = 0; end

if isempty(filtParam.fs), error('Sampling rate should be specified !');end
if isempty(filtParam.fmin)
    filtParam.fmin = 60.;
    warning('No low frequency limit. f = %4.2e used as default.',...
        paramFilt.fmin);
end
if isempty(filtParam.fmax)
    filtParam.fmax = 8000.;
    warning('No high frequency limit. f = %4.2e used as default.',...
        paramFilt.fmax);
end

if isempty(filtParam.noct)
    filtParam.noct = 1;
    warning('No octave band factor defined. noct = %d used as default.',...
        paramFilt.noct);
end

if filtParam.plotfl
    hh = figure('Name','Filter bank response');
    hold on;
    Hsum = 0;
    nbfreq = 2^15;
end


% center and cutoff frequencies of each filter
freq = freqoct(filtParam.noct,filtParam.fmin,filtParam.fmax,1);
% check Nyquist condition
idx_keep = freq(:,3) <= filtParam.fs/2;
freq = freq(idx_keep,:);
nbfilt = size(freq,1);

for ifilt = 1 : nbfilt
    filterBank(ifilt).fc = freq(ifilt,2);
    Ws = [freq(ifilt,1) freq(ifilt,3)];
    Wn = Ws / (filtParam.fs/2);
    sos = getFilterCoef(1,filtParam.order,Wn, filtParam.energyfl);
    
    filterBank(ifilt).sos = sos; 
    
    if filtParam.plotfl               
        [H,ff] = freqz(sos,nbfreq,filtParam.fs);
        Hsum = Hsum + H;
        Hdb = 20*log10(abs(H));
        
        plot(ff,Hdb,'linewidth',1);                
    end
    
    
end

if filtParam.plotfl
    Hsum_db = 20*log10(abs(Hsum));
    plot(ff,Hsum_db,'linestyle','-','linewidth',2);
    grid on
    hold off
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]');
    set(gca,'xlim',[freq(1,2)*2^(-1/filtParam.noct) filtParam.fs/2]);
    set(gca,'xScale','log','xTick',freq(:,2));
    set(gca,'ylim',[-80 10])
end



end

function [sos] = getFilterCoef(type,order,Wn,energyfl,rp,rs)

if nargin < 6, rs = []; end
if nargin < 5, rp = []; end

if type > 4, type = 1; end

switch(type)
    % Butterworth
    case 1
        [z,p,k] = butter(order,Wn,'bandpass');
        sos = zp2sos(z,p,k);
    case 2
        
    case 3
        
    case 4
        
end

if energyfl
    sos = sos.*sos;
end

end