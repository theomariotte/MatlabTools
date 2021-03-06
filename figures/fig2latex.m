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

nbl = length(ax);

colors = {[.2 .6 .8],[.8 .4 .3],'k','k','c'};
fprintf('Change figure format (font size, dimensions,...)\n');

% get axes dimensions (useful to resize axes object when deleting a sublot
% axes)
dummyfig = figure;
plot(1,1);
ax_pos = get(gca,'position');
delete(dummyfig);

for k = 1 : nbl
    % get lines on each axes
    hls = get(ax(k),'Children');
    
    if isempty(hls) || k == 1
        delete(ax(k))
    else            
        %%% Change axes properties
        % font size
        set(ax(k),'position',ax_pos);
        set(ax(k),'fontsize',14);
        % tick interpreter to latex
        set(ax(k),'TickLabelInterpreter','latex');
        set(ax(k),'Title',[]);
        %set(ax(k),'ylim',[-60 0])
        % tick labeling (automatic)
        %set(ax,'xTickLabelMode','auto');
        % ylimits
        %set(ax,'ylim',[-100 0]);
        % labels
        xlabel('Fr�quence [s]');
        ylabel('DSP [dBA]')        

        %%% change line properties
        % Line width
        set(hls,'linewidth',2);
        %set(hls(1),'Visible','off');
        set(hls(1),'color',colors{1});
        set(hls(2),'linestyle','--')
        set(hls(2),'color',colors{2});
        % line color
        %    set(hls,'color','k');

    end
end


<<<<<<< HEAD
=======
%%% change line properties
% Line width
set(hls,'linewidth',2);
% line color
% set(hls,'color','k');
>>>>>>> LatexFigDev

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