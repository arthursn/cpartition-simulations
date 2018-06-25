#!/bin/bash

# python3 plot_cprofiles_sep.py coupled_FoFo_375_CCE -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -special
# python3 plot_cprofiles_sep.py coupled_FoFo_375_CCEpara -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -special
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu30e3 -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -special

# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu16e3 -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -ylim -.1:2.1
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu20e3 -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -ylim -.1:2.1
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu23e3 -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -ylim -.1:2.1
# python3 plot_cprofiles_sep.py coupled_FoFo_375_CCEortho -show -t .1 1 10 100 1000 -xlim -1.16:1.16 -wbs -all -save -ylim -.1:2.1

# python3 plot_cprofiles.py mart_FoFo_CCE*py -t .1 1 10 100 1000 -xlim -1.16:0 -ylim -.1: -label -figsize 4.5 4.5 -save
python3 plot_cprofiles.py mart_FoFo_mu*py -t .1 1 10 100 1000 10000 -xlim -1.16:0 -ylim -.1: -label -figsize 4.5 4.5 -save -show

# python3 plot_cprofiles.py coupled_FoFo_375_CCE -t 1e-1 5.5 1e3 -label -show -tracking  -mirror -loc ll -save -suffix _tracking -xlim -.45:.1 -ylim -.1:1.9 -show