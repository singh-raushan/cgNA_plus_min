function [] = plot2D_cgDNAp_coordinates(w_groundstate, w_initial_guess, w_circle, do_title, seq_name, guess_nb, do_save, is_diff, varargin)
% Add Utilities to path
addpath('./cgDNAp_periodic');

% Check varargin for optional input
p = CheckInput(varargin);

nbp = max([length(w_groundstate),length(w_initial_guess),length(w_circle)])/24 ;

% Initialise the 3d array that will contains the coordinates
allCoordcgDNA = nan(nbp,12,3) ;
allCoordPho = nan(nbp,12,3) ;

% For loop over the number of sequence
for i  = 1:3
    % Read ground state
    switch i
        case 1
            shape = w_groundstate;
        case 2
            shape = w_initial_guess;
        case 3
            shape = w_circle;
    end

    % By default transfrom the groundstate in dimensional one
    if p.RotDeg
        shape = nondim2curPeriodic(shape);
    end

    % Decompose the shape into intra/inter and W/C pho rotations and
    % translations
    [Buckle_Propeller_Opening, ...
        Shear_Stretch_Stagger, ...
        pho_W_rot , pho_W_tr, ...
        Tilt_Roll_Twist, ...
        Shift_Slide_Rise, ...
        pho_C_rot , pho_C_tr ] = vector2shapesPeriodic(shape);

    % Collects intras
    intra = [Buckle_Propeller_Opening, ...
        Shear_Stretch_Stagger];
    % Collects inters
    inter = [Tilt_Roll_Twist, ...
        Shift_Slide_Rise];
    % Collects Watson Phosphate
    pho_W = [pho_W_rot, ...
        pho_W_tr ] ;
    % Collects Crick Phosphate
    pho_C = [pho_C_rot, ...
        pho_C_tr ] ;

    % Store the decomposed cgDNA coordinates into 3 dimensional array
    tmp = [ inter intra ];
    allCoordcgDNA(1:size(tmp,1),1:size(tmp,2),i) = tmp;

    % Store the decomposed Pho coordinates into 3 dimensional array
    tmp = [  pho_W pho_C ];
    allCoordPho(1:size(tmp,1),1:size(tmp,2),i)=tmp;

end

% Initialise the figure
hsize = get(0,'ScreenSize');
fig1=figure('Position',hsize);
% fig1 = figure();

% Define x axis limits and ticks
xlim = [0 nbp+1];
xtick = setdiff(sort([ 1 0:10:nbp]),0);
x = 1:nbp;

% For loop over the 12 panels for cgDNA coordinates
for j=1:12

    subplot(4,3,j);

    % Read the data to plot in jth panel
    coor = squeeze(allCoordcgDNA(:,j,:));

    % Plot the data
    plot(x, coor,'LineWidth', p.LineWidth,'Marker',p.Marker,'MarkerSize',p.MarkerSize);

    % Add title, x, and y label
    set(gca,'XLim',xlim,'XTick',xtick,... %'xticklabels',{'A','T','A','T','A','T','A','T','A','T','A','T/A','A','T','A','T','A','T',},
        'Fontsize',p.FontSize);
    ylabel(p.ylabel{j}) ;
    xlabel(p.xlabel{j});
    title(p.title{j});

    grid on
end

% Show the legend
if not(is_diff)
    lgtxt =  {'Periodic groundstate (linear)', 'Initial guess (minicircle) ', 'Minicircle (final min energy shape)'};
    %%%% lgtxt =  {'Ground state', 'Initial guess bBDNA', 'Mini-circle'};
    legend1=legend(lgtxt) ;
    set(legend1,...
        'Position', [0.584491221110026,0.935386503530833,0.31953125,0.039968279143537],... [0.462211608886719 0.950662739322533 0.0901321411132813 0.0360824742268041],...
        'FontSize',p.FontSize,...
        'NumColumns',3);
end


% Set up axes position
axes('Position',[0,0,1,1],'visible','off');
% set title
if not(is_diff)
    if do_title
        fig1_title = "Inter and Intra Coordinates";%, sequence "+seq_name;
        %  fig1_title = "cgDNA coordinates" + seq_name+ ", link " + guess_nb;
        % fig1_title = seq_name + ", link " + guess_nb + ": difference in cgDNA coordinates";
    else
        fig1_title = "Difference in cgDNA coordinates";
    end
else
    fig1_title = "Difference: groundstate - circle shape, cgDNA coordinates";
end
sgtitle(fig1_title,'EdgeColor','black');

hold off

% Link all the panels for allowing simulatanous zooming
all_ha = findobj( fig1, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );

% save figure
if do_save
    savename1 = "../Results/2D_plots/"+seq_name+"_2Dsol"+guess_nb+"_PCcgDNA+min_cgDNA-coord";
    %    savefig(fig1, savename1);
end

% Initialise the figure --------------------
hsize = get(0,'ScreenSize');
fig2=figure('Position',hsize);
% fig2=figure();

% Define x axis limits and ticks
xlim = [0 nbp+1];
xtick = setdiff(sort([ 1 0:10:nbp]),0);
x = 1:nbp;

% For loop over the 12 panels for Pho coordinates
for j=1:12

    subplot(4,3,j);

    % Read the data to plot in jth panel
    coor = squeeze(allCoordPho(:,j,:));

    % Plot the data
    plot(x, coor,'LineWidth', p.LineWidth,'Marker',p.Marker,'MarkerSize',p.MarkerSize);

    % Add title, x, and y label
    set(gca,'XLim',xlim,'XTick',xtick,... %,'xticklabels',{'A','T','A','T','A','T','A','T','A','T','A','T/A','A','T','A','T','A','T',}
        'Fontsize',p.FontSize);
    ylabel(p.ylabel{j}) ;
    xlabel(p.xlabel_pho{j});
    title(p.title_pho{j});

    grid on
end

%%%% Show the legend
if not(is_diff)
    lgtxt =  {'Periodic groundstate (linear)', 'Initial guess (minicircle) ', 'Minicircle (final shape)'};
    %%%% lgtxt =  {'Periodic groundstate (linear)', 'Initial guess (minicircle) ', 'Mini-circle'};
    legend2=legend(lgtxt) ;
    set(legend2,...
        'Position', [0.584491221110026,0.935386503530833,0.31953125,0.039968279143537],... [0.462211608886719 0.950662739322533 0.0901321411132813 0.0360824742268041],...
        'FontSize',p.FontSize,...
        'NumColumns',3);
end

% Set up axes position
axes('Position',[0,0,1,1],'visible','off');
% set title
if not(is_diff)
    if do_title
        fig2_title = "Phosphate Coordinates"; %, sequence "+seq_name; %+", guess " + guess_nb;
        % fig2_title = seq_name + ", link " + guess_nb + ": difference in Phosphate coordinates";
    else
        fig2_title = "Difference in Phosphate coordinates";
    end
else
    fig2_title = "Difference: groundstate - circle shape, Phosphate coordinates";
end
sgtitle(fig2_title,'EdgeColor','black');

hold off

% Link all the panels for allowing simulatanous zooming
all_ha = findobj( fig2, 'type', 'axes', 'tag', '' );
linkaxes( all_ha, 'x' );

% save figure
if do_save
    savename2 = "../Results/2D_plots/"+seq_name+"_2Dsol"+guess_nb+"_PCcgDNA+min_Pho-coord";
    %    savefig(fig2, savename2);
end

end

%--------------------------------------------------------------------------

function p = CheckInput(inputarg)

p = inputParser ;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validMarker = @(x) ischar(x) && contains('.ox+*sdv^<>phn',x);
validUnit   = @(x) isnumeric(x) && x == 1 || x == 0;

addOptional(p,'RotDeg',1,validUnit) ;
addOptional(p,'FontSize',15,validScalarPosNum) ;
addOptional(p,'LineWidth',1,validScalarPosNum) ;
addOptional(p,'Marker','o',validMarker) ;
addOptional(p,'MarkerSize',1,validMarker) ;

parse(p,inputarg{:});

p=p.Results;

if p.RotDeg

    p.ylabel = {
        'degrees'; 'degrees'; 'degrees'; ...
        'angstrom'; 'angstrom'; 'angstrom'; ...
        'degrees'; 'degrees'; 'degrees'; ...
        'angstrom'; 'angstrom'; 'angstrom' ...
        };
else
    p.ylabel = {
        'rad/5'; 'rad/5'; 'rad/5'; ...
        'angstrom'; 'angstrom'; 'angstrom' ; ...
        'rad/5'; 'rad/5'; 'rad/5'; ...
        'angstrom'; 'angstrom'; 'angstrom' ...
        };
end

p.title = {
    'Tilt'; 'Roll'; 'Twist'; ...
    'Shift'; 'Slide'; 'Rise'
    'Buckle'; 'Propeller'; 'Opening'; ...
    'Shear'; 'Stretch'; 'Stagger' ...
    };

p.title_pho = {
    'W Rot 1'; 'W Rot 2'; 'W Rot 3'; ...
    'W Tra 1'; 'W Tra 2'; 'W Tra 3' ; ...
    'C Rot 1'; 'C Rot 2'; 'C Rot 3'; ...
    'C Tra 1'; 'C Tra 2'; 'C Tra 3' ...
    };

p.xlabel = { 'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair'; 'basepair'; 'basepair'; ...
    'basepair'; 'basepair'; 'basepair' ...
    };

p.xlabel_pho = { 'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'; ...
    'basepair step'; 'basepair step'; 'basepair step'...
    };

end