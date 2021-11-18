
ROOT_URL='https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/include'

for F in \
angle-utils.hpp \
angle.hpp \
character.hpp \
ellipsoid-wgs84.hpp \
ellipsoid.hpp \
math-const.hpp \
number.hpp \
sin_cos.hpp
do
	if [[ ! -f "$F" ]]
	then
		echo "$F"
		curl -O "$ROOT_URL/$F"
	fi
done
