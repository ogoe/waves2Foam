function processData( toolpath )

close all, clear all, clc

homeDir = pwd;

if nargin == 0
    addpath(sprintf('%s/../../../../applications/utilities/misc/matlab/postprocessing',pwd));
else
    addpath( toolpath );
end

%% Load the numerical simulation
cd ../surfaceElevation/
h = dir('0*');
cd ..

[time, x, y, z, eta] = readSurfaceElevation( h(1).name );

cd ( homeDir );


%% Load the experimental data
fileNames  = {'a04.dat';'a10.dat';'a13.dat';'a14.dat';'a15.dat';'a17.dat';'a19.dat';'a21.dat'};
xCoord     = [4.0; 10.5; 13.5; 14.5; 15.7; 17.3; 19.0; 21.0];

expData    = cell(length(fileNames),1);

for i=1:length(expData)
    expData{i} = load(sprintf('../experimentalData/%s',fileNames{i}));
end

%% Plot data

tOffset = 0.53 - 2.02;

hf = figure; figureSize(1.4,2.5)

plotLegend = true;
fontsize = 17;

for i=1:length(expData)
    subplot(2,4,i), hold on, grid on
    title(sprintf('x = %.1f m',xCoord(i)),'fontname','times','fontsize',fontsize);
    
    I = find(x == xCoord(i));
    
    if ~isempty(I) && length(I) == 1
        plot(time - tOffset, eta(:,I),'linewidth',1.)
        plot(expData{i}(:,1), expData{i}(:,2),'--r','linewidth',1.)
        
        
        if plotLegend
            legend('Simulation','Experiment','Location','NorthWest')           
            
            plotLegend = false; 
        end
        
        
    end
    
    set(gca,'fontname','times','fontsize',fontsize,'box','on','xlim',[30 40],'ylim',[-0.02 0.03]);

    if i==1 || i == 5
        ylabel('Elevation, [m]');
    end
    
    if i > 4
        xlabel('Time, [s]');
    end
end

fprintf('\n\n- Printing the figure to <case>/matlab/result.eps\n\n');
print -depsc result

time
eta

function figureSize(h,w)
% Rescaling the plotting window

hpp = get(gcf,'PaperPosition');

hpp(3) = w * hpp(3);
hpp(4) = h * hpp(4);

set(gcf,'PaperPosition',hpp);

hp = get(gcf,'Position');

hp(3) = w * hp(3);
hp(4) = h * hp(4);

set(gcf,'Position',hp);