close all, clear all

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

xx = load('../postProcessedWaves2Foam/surfaceElevation/surfaceElevation_indexXYZ.txt');

for i=1:size(xx,1)
    data = load(sprintf('../postProcessedWaves2Foam/surfaceElevation/surfaceElevation_%d.txt',xx(i)));
    
    if i==1
        eta = zeros(size(xx, 1), size(data, 1));
    end
    
    eta(i,:) = data';
end

%% Plot the comparison with experimental data
fSize = 14;

for i=0:2:22
    x = load(sprintf('../experimentalData/data_%02d.dat',i));
    
    subplot(6,2,i/2+1), hold on, grid on, set(gca,'box','on');
    fill([30 59 59 30]/100,[0 0 27 27]/100,0.75*[1 1 1])
    
    plot(xx(:,2), eta(:,i/2+1), 'linewidth',1.5)
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