#!/bin/bash



percentBar ()  { 
    local prct totlen=$((8*$2)) lastchar barstring blankstring;
    printf -v prct %.2f "$1"
    ((prct=10#${prct/.}*totlen/10000, prct%8)) &&
        printf -v lastchar '\\U258%X' $(( 16 - prct%8 )) ||
            lastchar=''
    printf -v barstring '%*s' $((prct/8)) ''
    printf -v barstring '%b' "${barstring// /\\U2588}$lastchar"
    printf -v blankstring '%*s' $(((totlen-prct)/8)) ''
    printf -v "$3" '%s%s' "$barstring" "$blankstring"
}


dir=$1

foo=(`ls $dir/*.fasta`)

for index in "${!foo[@]}"; do

    fullfilename=${foo[$index]}
    len=${#foo[@]}

    i=$(( 100 * index / len ))

    filename=$(basename "$fullfilename")
    fname="${filename%.*}"
    ext="${filename##*.}"

#     echo "Dir: $dir"
#     echo "Input File: $fullfilename"
#     echo "Filename without Path: $filename"
#     echo "Filename without Extension: $fname"
#     echo "File Extension without Name: $ext"
 
    awk -v var="${fname%%_*}" -i inplace '/^>/{print ">" ++i "_" var; next}{print}' "$fullfilename"


#     TODO: Add a beautiful progress bar
    percentBar $i $((COLUMNS-50)) bar
    printf '\r\e[40;44m%s\e[0m%6.2f%%' "$bar" $i
    read -srt .002 _ && break
done
cat $dir/*.fasta > "$dir/combined.fasta"