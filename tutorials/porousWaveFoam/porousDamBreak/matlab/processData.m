close all, clear all

addpath('./../../../../applications/utilities/misc/matlab/postprocessing/');

figure
% Change the size of the figure
w = 1; h = 2
hpp = get(gcf,'PaperPosition');
hpp(3) = w * hpp(3);
hpp(4) = h * hpp(4);
set(gcf,'PaperPosition',hpp);

hp = get(gcf,'Position');
hp(3) = w * hp(3);
hp(4) = h * hp(4);
set(gcf,'Position',hp);

%% Read simulation data

cd ../postProcessing
[time, xx, y, z, eta] = readSurfaceElevation('0');
cd ../matlab

%% Plot the comparison with experimental data
fSize = 14;

for i=0:2:22
    x = load(sprintf('../experimentalData/data_%02d.dat',i));
    
    subplot(6,2,i/2+1), hold on, grid on, set(gca,'box','on');
    fill([30 59 59 30]/100,[0 0 27 27]/100,0.75*[1 1 1])
    
    dt = time - (i/2)*0.2;
    It = find(abs(dt) == min(abs(dt)));
    plot(xx, eta(It,:), 'linewidth',1.5)
    plot(x(:,1), x(:,2), 'bo')
    
    xlim([0 0.892]), ylim([0 0.3])
    
    if mod(i/2+1,2)==0
        set(gca,'yticklabel',{''});
        pos = get(gca,'position');
        set(gca,'position', pos - [0.05 0 0 0]);
    else
        ylabel('y, [m]','fontsize',14);
    end
    
    if i < 20
        set(gca,'xticklabel',{''});
    else
        xlabel('x, [m]','fontsize',14);
    end
    
    set(gca,'fontsize', fSize)
    text(0.02, 0.05,sprintf('t = %.1f s',i/10), 'fontsize', fSize - 2);
end

print -depsc comparison