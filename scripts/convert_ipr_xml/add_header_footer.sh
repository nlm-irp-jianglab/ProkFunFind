a=$1
echo -e "$(cat header)\n$(cat $a)\n$(cat footer)" >${a}
