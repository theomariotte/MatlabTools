function [] = fig2latex(figname,outFmt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure formating for LaTeX from .fig file
% formated figure will be saved in the same directory
% Inputs :
%   - figname : path to the figure (including name and extension)
%   - outFmt : extension of the generated picture
% Available out format :
% 'fig', 'eps','png','pdf'
%
% WARNING : for the moment, this function does not take multiple axes plots
% into account (e.g. figures generated using subplot)
% @Author T. Mariotte
% @todo (25/05/2020) add parameters structure in input argument to define
% parameters such as scales, axes limits and so on ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2, outFmt = []; end
if isempty(outFmt)
    outFmt = 'eps';
    warning('Default output format used : %s',outFmt);
end

[figpath, fname, ext] = fileparts(figname);
if strcmp(ext,'.fig')==0, error('Input file should be ''.fig'' format'); end

% open figure
hfig = openfig(figname);

% get axes of the figure
ax = hfig.Children;

% get lines on each axes
hls = get(ax,'Children');

fprintf('Change figure format (font size, dimensions,...)\n');

%%% Change axes properties
% font size
set(ax,'fontsize',14);
% tick interpreter to latex
set(ax,'TickLabelInterpreter','latex');
% tick labeling (automatic)
set(ax,'xTickLabelMode','auto');
% ylimits
set(ax,'ylim',[-100 0]);
% labels
xlabel('Fréquence [Hz]');
ylabel('Amplitude [dB]')

%%% change line properties
% Line width
set(hls,'linewidth',2);
% line color
% set(hls,'color','k');

%%% Output formating ('eps' recommanded for LaTeX use)
switch outFmt
    case 'eps'
        fmt = '-depsc';
    case 'png'
        fmt = '-dpng';
    case 'jpg'
        fmt = '-djpeg';
    case 'pdf'
        fmt = '-dpdf';
    otherwise 
        error('No file format found. Try : ''eps'', ''png'',''jpg'',''pdf''');
end

%%% Save the new figure i, the same directory
print([figpath fname],fmt);
fprintf('New figure saved as :\n %s\n',[figpath fname '.' outFmt]);


end