#!/usr/bin/bash

if [[ ! -f angle.hpp ]]
then
	curl -O 'https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/include/angle.hpp'
fi

if [[ ! -f angle-utils.hpp ]]
then
	curl -O 'https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/include/angle-utils.hpp'
fi

if [[ ! -f ellipsoid.hpp ]]
then
	curl -O 'https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/include/ellipsoid.hpp'
fi

if [[ ! -f ellipsoid-wgs84.hpp ]]
then
	curl -O 'https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/include/ellipsoid-wgs84.hpp'
fi
