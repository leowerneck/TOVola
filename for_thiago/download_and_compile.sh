#!/bin/sh

# Error out in case of a problem
set -e

num_cores=16
tov_th_url="https://raw.githubusercontent.com/leowerneck/TOVola/refs/heads/shift/for_thiago/tov.th"
tov_par_url="https://raw.githubusercontent.com/leowerneck/TOVola/refs/heads/shift/for_thiago/tov.par"

# Download GetComponents using instructions from einsteintoolkit.org
curl -kLO "https://raw.githubusercontent.com/gridaphobe/CRL/ET_2025_05/GetComponents"
chmod a+x GetComponents

# Download tov.th from Leo's TOVola repo
curl -kLO "$tov_th_url"

# Download the toolkit
./GetComponents --shallow tov.th --root Cactus_tov

# Download the TOV parfile
cd Cactus_tov
curl -kLO "$tov_par_url"

# Compile the toolkit
./simfactory/bin/sim setup-silent
./simfactory/bin/sim build -j$num_cores --thornlist=thornlists/tov.th

# Cleanup
rm -f GetComponents
