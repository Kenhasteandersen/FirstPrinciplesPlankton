Code to produce the figures in "From first principles to ecosystem structure andfunction of unicellular plankton":

R code to produce the figures:

 * makefigures.R. Code to do almost all figures. run "plotAll()".
 * app.R. Shiny app with the chemostat setup.
 * modelChemostat.R. Code for the chemostat runs
 * model.R. Core routines for the unicellular model.
 * basetools.R. Plotting routines.
 * matlab: code to run the two figures with the water column and global setups
 * lib: Pre-compiled versions of the NUMmodel library needed to do the simulations. If they do not work then get the source code https://github.com/Kenhasteandersen/NUMmodel/releases/tag/v0.91 and compile yourself.