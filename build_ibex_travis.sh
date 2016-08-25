#!/bin/bash

set -x
sudo chmod a+rw /opt/
# check to see if ibex folder is empty
if [ ! -e "/opt/ibex/lib/libibex.a" ]; then
	cd /opt/
	git clone https://github.com/ibex-team/ibex-lib.git
	cd ibex-lib
	./waf configure
	sudo ./waf install -j2
	export PKG_CONFIG_PATH=/usr/local/share/pkgconfig/
else
  echo 'Using cached directory.';
fi

# check to see if ibex-geometry folder is empty
if [ ! -e "/opt/ibex-geometry/LICENSE" ]; then
	cd /opt/
	git clone https://github.com/benEnsta/ibex-geometry.git
	cd ibex-geometry
	cmake .
	sudo make install -j2
else
  echo 'Using cached directory.';
fi

# check to see if ibex-geometry folder is empty
if [ ! -e "/opt/VIBES/LICENSE" ]; then
	cd /opt/
	git clone https://github.com/ENSTABretagneRobotics/VIBES.git
else
  echo 'Using cached directory.';
fi

