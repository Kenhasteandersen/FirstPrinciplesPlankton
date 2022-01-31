<<<<<<< HEAD
addpath('~/NUMmodel/matlab')

p = setupGeneralistsOnly(20,true);
p = parametersGlobal(p,2);
=======
%p = setupGeneralistsOnly(20,true);
%p = parametersGlobal(p,2);
>>>>>>> b4a02ff95bf283b443b333257724866e92e9bdcc

%sim = simulateGlobal(p);

load ecco20

%%
figure(1)
[tiles,sim] = plotGlobalBiomass(sim, 'eckert4');
if exist('setFigWidth')
    setFigWidth(10)
    setFigHeight(14)
    exportgraphics(tiles, '../globalbiomass.pdf')
end
%print -dpdf globalbiomass.pdf

figure(2)
[sim, tiles] = plotGlobalFunctions(sim,'eckert4');
if exist('setFigWidth')
    setFigWidth(10)
    setFigHeight(14)
    exportgraphics(tiles, '../globalfunctions.pdf')
end

fprintf("Total gross production %f.2 Pg_C/yr\n", mean(sim.ProdGrossTotal*1e-15));
fprintf("Total net production   %f.2 Pg_C/yr\n", mean(sim.ProdNetTotal*1e-15));
fprintf("Total htl production   %f.2 Pg_C/yr\n", mean(sim.ProdHTLTotal*1e-15));

%print -dpdf globalfunctions.pdf

%save ecco20 sim
