%
% Make three water column runs
%
addpath('~/Documents/Source/NUMmodel/matlab')
addpath('~/Documents/Source/Matlab/')

p = setupGeneralistsOnly(50);
p = parametersWatercolumn(p);
p.tEnd = 4*365;
setHTL(0.1);
p.dt=0.01;

days = [120 200];
% figure(1)
% sim = simulateWatercolumn(p,0,330);
% plotOneWatercolumn(sim); % Eutrophic
% 
% figure(2)
% sim = simulateWatercolumn(p,20,306);
%plotOneWatercolumn(sim); % Oligotrophic

%%
sim = simulateWatercolumn(p, 60,-10);
p.tEnd = 365;
sim = simulateWatercolumn(p, 60,-10,sim);
%%
figure(3)
clf
t = tiledlayout(1,3,'TileSpacing','compact','padding','compact');
plotWatercolumnTime(sim,depthMax=200,bNewPlot=false)
nexttile(1)
title('')
plotlabel('a)',false)
%set(gca,'XTickLabel','')
xlabel('')
defaultAxes(16)

nexttile(2)
title('')
plotlabel('b)',false)
set(gca,'YTickLabel','')
ylabel('')
defaultAxes()

nexttile(3)
title('')
plotlabel('c)',false)
ylabel('')
set(gca,'YTickLabel','')
%caxis([0.1 300])
xlabel('')
defaultAxes()

for i = 1:3
    nexttile(i)
    line(days(1)*[1 1],[-200 0],1000*[1,1],'color','w','linestyle','--','linewidth',2)
    line(days(2)*[1 1],[-200 0],1000*[1,1],'color','w','linestyle','--','linewidth',2)
end

xlabel(t,'Time (days)')
setFigHeight(4.5)

%set(gca, 'xtick',1:30:365);%,...
    %'xticklabel',{'Jan','Apr','Jul','Oct'})

exportgraphics(t, '../watercolumn1.pdf')
%%
figure(4)
t = tiledlayout(2,2,'TileSpacing','compact','padding','compact',...
    'tileindexing','columnmajor');
%nexttile
plotWatercolumn(sim, days(1), depthMax=200, bNewplot=false)

%nexttile
plotWatercolumn(sim, days(2), depthMax=200, bNewplot=false)

nexttile(1)
title('')

colorbar off
nexttile(2)
set(gca,'yticklabel',-200:50:0)
nexttile(3)
title('')

slabel='dfeg';
for i = 1:4
    nexttile(i)
    
    ylabel('')
    xlabel('')
    if (i<3)
        %set(gca,'yticklabel',-200:50:0)
    else
        set(gca,'YTickLabel','')
    end
    set(gca,'ytick',-200:50:0)
    
    plotlabel(strcat(slabel(i),')'), false)
    defaultAxes(16)
end

xlabel(t,'Cell mass ({\mu}g_C)')
ylabel(t,'Depth (m)')

setFigHeight(10)

exportgraphics(t, '../watercolumn2.pdf')

function plotOneWatercolumn(sim)
    
    clf
    ylimit = [-200 0];
    
    ixTime = 3*365:4*365;
    tiledlayout(1,4,'TileSpacing','compact','padding','compact')
    
    nexttile
    N = squeeze(mean(sim.N(:,ixTime),2))/14; % Convert units to umolN/l
    plot(N,-sim.z,'b-','linewidth',2)
    hold on
    
    DOC = squeeze(mean(sim.DOC(:,ixTime),2))/12;% Convert units to pmolN/l
    plot(DOC,-sim.z,'color','#964B00','linewidth',2)
    
    %L = squeeze(mean(sim.L(:,ixTime),2));
    %plot(L, -sim.z, 'y-','linewidth',2)
    
    ylim(ylimit)
    legend({'DIN ({\mu}mol N/l)', 'DOC (pmol C/l)', 'Light ({\mu}E/m^2/day)'})
    
    plotWatercolumn(sim,365*3.5, depthMax=-min(ylimit), bNewPlot=false)

    nexttile
    sim.B = mean(sim.B(5,:,ixTime),3);
    panelSpectrum(sim, 1)
    
end
