function filterBank = designEqFilt(param)
%% Banc de filtres passe-bande - décomposition en 1/n d'octave
% 
% param : 
%   - 'type' : type of filter 
%              (1) Butterworth
%              (2) Elliptic
%              (3) Chebychev I
%              (4) Chebychev II
%   - 'rp' : pass band ripple (dB)
%   - 'rs' : stop band attenuation (dB)
%   - 'fmin' : lowest band center frequency
%   - 'fmax' : 'highest band center frequency
%   - 'fs'   : sampling rate
%   - 'noct' : band decomposition (e.g. noct = 3 means 1/3 octave bands)
%   - 'lateralBands' : include or not Low pass and High pass filters for
%   lowest and highest bands
%              (0) no lateral bands
%              (1) add them
%   - 'energyfl' : double filtering to conserve energy (useful for
%   equalization)
%   - 'order' : filter order
%   - plotfl : (0) no plot ; (1) plot frequency response
%
% Matlab tools/filter

if ~isfield(param,'type'), param.type = []; end
if ~isfield(param,'rs'), param.rs = 40.; end
if ~isfield(param,'rp'), param.rp = 1.; end
if ~isfield(param,'fmin'), param.fmin = []; end
if ~isfield(param,'fmax'), param.fmax = []; end
if ~isfield(param,'fs'), param.fmin = []; end
if ~isfield(param,'noct'), param.noct = []; end
if ~isfield(param,'order'), param.order = 6; end
if ~isfield(param,'lateralBands'), param.lateralBands = 0; end
if ~isfield(param,'energyfl'), param.energyfl = 0; end
if ~isfield(param,'plotfl'), param.plotfl = 0; end

if isempty(param.fs), error('Sampling rate should be specified !');end
if isempty(param.fmin)
    param.fmin = 60.;
    warning('No low frequency limit. f = %4.2e used as default.',...
        paramFilt.fmin);
end
if isempty(param.fmax)
    param.fmax = 8000.;
    warning('No high frequency limit. f = %4.2e used as default.',...
        paramFilt.fmax);
end

if isempty(param.noct)
    param.noct = 1;
    warning('No octave band factor defined. noct = %d used as default.',...
        paramFilt.noct);
end
if isempty(param.type)
    param.type = 1;
    warning('No filter type found. Default = Butterworth');
end

if param.plotfl
    hh = figure('Name','Filter bank response');
    hold on;
    Hsum = 0;
    nbfreq = 2^15;
end

% In case of LP/HP filters for lateral bands
lateral = param.lateralBands;
if lateral
    param.fmin = param.fmin *2^(-1/param.noct);
    param.fmax = param.fmax *2^(1/param.noct);
end

% center and cutoff frequencies of each filter
freq = freqoct(param.noct,param.fmin,param.fmax,1);
% check Nyquist condition
idx_keep = freq(:,3) <= param.fs/2;
freq = freq(idx_keep,:);
nbfilt = size(freq,1);

for ifilt = 1 : nbfilt
    
    if lateral==1 && ifilt==1
        cut_type = 'low';
        Ws = freq(ifilt,3);
    elseif lateral == 1 && ifilt == nbfilt        
        cut_type = 'high';
        Ws = freq(ifilt,1);
    else
        cut_type = 'bandpass';
        Ws = [freq(ifilt,1) freq(ifilt,3)];
    end
    
    filterBank(ifilt).fc = freq(ifilt,2);
    Wn = Ws / (param.fs/2);
    sos = getFilterCoef(param.type,param.order,Wn,...
        param.rp,...
        param.rs,...
        cut_type);
    
    filterBank(ifilt).sos = sos; 
    filterBank(ifilt).energyfl = param.energyfl;
    
    if param.plotfl               
        [H,ff] = freqz(sos,nbfreq,param.fs);
        if param.energyfl
            H = H.*H;
        end
        Hsum = Hsum + H;
        Hdb = 20*log10(abs(H));
        
        plot(ff,Hdb,'linewidth',1);                
    end
    
    
end

if param.plotfl
    Hsum_db = 20*log10(abs(Hsum));
    plot(ff,Hsum_db,'linestyle','-','linewidth',2);
    grid on
    hold off
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]');
    set(gca,'xlim',[freq(1,2)*2^(-1/param.noct) param.fs/2]);
    set(gca,'xScale','log','xTick',freq(:,2));
    set(gca,'ylim',[-80 10])
end



end

function [sos] = getFilterCoef(type,order,Wn,rp,rs,cut)

if type > 4, type = 1; end

switch(type)
    % Butterworth
    case 1
        if isempty(order)
            error('No order found !');
            % TODO : get order from frequencies
            [order,Wn] = buttord(Wn,Wn,rp,rs);
        end
        [z,p,k] = butter(order,Wn,cut);  
        
    % elliptic
    case 2
        if isempty(order)
            error('No order found !');
            [order,Wn] = ellipord(Wn(1),Wn(2),rp,rs);
        end
        [z,p,k] = ellip(order,rp,rs,Wn,cut);
        
    % Cheby I
    case 3
        if isempty(order)
            error('No order found !');
            [order,Wn] = cheb1ord(Wn(1),Wn(2),rp,rs);
        end
        [z,p,k] = cheby1(order,rp,Wn,cut);
        
    % Cheby II
    case 4
        if isempty(order)
            error('No order found !');
            [order,Wn] = cheb2ord(Wn(1),Wn(2),rp,rs);
        end
        [z,p,k] = cheby2(order,rp,Wn,cut);
end
sos = zp2sos(z,p,k);


end