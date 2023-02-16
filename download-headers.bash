
# header files
for F in \
ellipsoid-wgs84.hpp \
ellipsoid.hpp
do
	if [[ ! -f "$F" ]]
	then
		echo "$F"
		curl -O https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/include/"$F"
	fi
done
