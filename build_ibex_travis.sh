#!/bin/bash

set -x
# check to see if ibex folder is empty
if [ ! -e "$HOME/ibex/lib/libibex.a" ]; then
	git clone https://github.com/ibex-team/ibex-lib.git
	cd ibex-lib
	./waf configure
	sudo ./waf install -j2
	export PKG_CONFIG_PATH=/usr/local/share/pkgconfig/
else
  echo 'Using cached directory.';
fi

# check to see if ibex-geometry folder is empty
if [ ! -e "$HOME/ibex-geometry/LICENSE" ]; then
	git clone https://github.com/benEnsta/ibex-geometry.git
	cd ibex-geometry
	cmake .
	sudo make install -j2
else
  echo 'Using cached directory.';
fi


