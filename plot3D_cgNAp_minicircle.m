function [] = plot3D_cgNAp_minicircle(shape, sequence, do_title, seq_name, guess_nb, do_save, varargin)
%--------------------------------------------------------------------------

% Check varargin for optional input
p = CheckInput(varargin);

% Add the folder containing the Viewer to the path
addpath(genpath('./cgDNApviewer'));
addpath('./cgDNAp_periodic');


% Initialise the figure
fig = figure();

cgDNApviewerPeriodic(shape, sequence, p.Reframe) ;

% By default do not plot axis
axis(p.flag_axis)

% figure title
if do_title
    fig_title = "cgNA+min minicircle"; %, sequence "+seq_name; %+", guess " + guess_nb;
    title(fig_title)
end

hold off

% save figure
if do_save
    savename = "../Results/3D_views/"+seq_name+"_3Dsol"+guess_nb+"_PCcgDNA+min";
    %    savefig(fig, savename)
end

end

function p = CheckInput(inputarg)

p = inputParser ;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validUnit   = @(x) isnumeric(x) && x == 1 || x == 0;

addOptional(p,'Reframe',1,validScalarPosNum) ;
addOptional(p,'Axis',0,validUnit) ;


parse(p,inputarg{:});

p=p.Results;

if p.Axis == 0
    p.flag_axis = 'off';

elseif p.Axis == 1
    p.flag_axis = 'on';

end

end

