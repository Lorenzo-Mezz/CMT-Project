/usr/local/bin/matlab-2021b -batch "run('src/tabletempandco2.m')"
gcc ./src/projet.c -o ./bin/projet -lm && ./bin/projet
/usr/local/bin/matlab-2021b -batch "run('src/plot_graphs.m')"
