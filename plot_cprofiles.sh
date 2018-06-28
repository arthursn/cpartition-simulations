#!/bin/bash

# Plot C profiles separately
# python3 plot_cprofiles_sep.py coupled_FoFo_375_CCE -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --wbs --all --special
# python3 plot_cprofiles_sep.py coupled_FoFo_375_CCEpara -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --wbs --all --special
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu30e3 -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --wbs --all --special
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu16e3 -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --ylim -.1 2.1 --wbs --all
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu20e3 -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --ylim -.1 2.1 --wbs --all
# python3 plot_cprofiles_sep.py coupled_FoFo_375_mu23e3 -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --ylim -.1 2.1 --wbs --all
# python3 plot_cprofiles_sep.py coupled_FoFo_375_CCEortho -t 1e-1 1e0 1e1 1e2 1e3 --xlim -1.16 1.16 --ylim -.1 2.1 --wbs --all

# Plot C profiles in a single figure
# python3 plot_cprofiles.py mart_FoFo_mu*py mart_FoFo_CCE*py -t 1e-1 1e0 1e1 1e2 1e3 1e4 --xlim -1.16 0 --label --figsize 4.5 4.5 --save
# python3 plot_cprofiles.py coupled_FoFo_375_CCE -t 1e-1 5.5 1e3 --xlim -.45 .1 --ylim -.1 1.9 --label --tracking --mirror --loc ll --append _tracking --save

# Plot muC profiles
python3 plot_muprofiles.py coupled_FoFo_375_CCE -t 1e-1 1e0 1e1 1e3 1e3 --xlim -1.16 1.16 --label --tracking --mirror --loc ll --append _tracking
# python3 plot_muprofiles.py coupled_FoFo_375_CCEortho -t 1e-1 1e0 1e1 1e3 1e3 --mu 8189 --xlim -1.16 1.16 --label --tracking --mirror --loc ll --append _tracking
# python3 plot_muprofiles.py coupled_FoFo_375_mu30e3 -t 1e-1 1e0 1e1 1e3 1e3 --mu 30e3 --xlim -1.16 1.16 --label --tracking --mirror --loc ll --append _tracking