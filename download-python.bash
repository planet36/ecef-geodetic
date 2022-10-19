
# scripts
for F in \
Nd-arange \
plot-points.py
do
	if [[ ! -f "$F" ]]
	then
		echo "$F"
		curl -O https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/bin/"$F"
	fi
done

# modules
for F in \
arange.py \
grouper.py \
remove_exponent.py
do
	if [[ ! -f "$F" ]]
	then
		echo "$F"
		curl -O https://raw.githubusercontent.com/planet36/dotfiles/main/link/.local/lib/python/"$F"
	fi
done
