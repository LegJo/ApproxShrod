#!/bin/bash

borneInf=$1
borneSup=$2
a=$3
b=$4
TypePot=$5
V0=$6
V1=$7
																		#fonction de traçage de graphique avec gnuplot
TracGraphGnuplot() {  													#1er arg = nom fichier, 2eme arg=ylabel, 3eme arg=xlabel, 4eme =min(initial $1), 5eme= max(initial $2)
    max=$(cut -d' ' -f2 "$1.tmp" | awk '{if($1>a)a=$1;}END{print a}')
    echo 'set term jpeg size 800,600' > gnu-script
    echo 'set offset 0,0,1,0' >> gnu-script
    echo 'set key top left' >> gnu-script
    echo 'set output "FctOnde.jpeg"' >> gnu-script
    echo 'set xrange ['$6':'$7']' >> gnu-script
    echo 'set xlabel "'$4'"' >> gnu-script
    echo 'set yrange [*:*]' >> gnu-script
    echo 'set ylabel "'$5'"' >> gnu-script
    echo "set title \"Fonction D'Onde z0="$8"  E="$9"\""  >> gnu-script 

    case $TypePot in
        1 ) 
            echo 'set arrow from -'$a',0 to -'$a','$max' nohead linecolor "red"' >> gnu-script
            echo 'f(x) = (-'$a'<=x && x<='$a') ? 0 : NaN' >> gnu-script
            echo 'set arrow from '$a',0 to '$a','$max' nohead linecolor "red"' >> gnu-script
            ;;
        2 )	
            echo 'f(x) = ('$6'<=x && x<=-'$b') ? 0 : (-'$b'<=x && x<='$b') ? '$V0' : ('$b'<=x && x<='$7') ? 0 : NaN' >> gnu-script
            ;;
        3 )
            echo 'f(x) = (-'$a'<=x && x<='$b') ? 0 : ('$b'<=x && x<='$a') ? '$V0' : NaN' >> gnu-script
            ;;
        4 ) 
            echo 'set arrow from -'$a',0 to -'$a','$max' nohead linecolor "red"' >> gnu-script
            echo 'f(x) = (-'$a'<=x && x<=-'$b') ? 0 : (-'$b'<=x && x<='$b') ? (-'$V0' * x**2 + '$V1') : ('$b'<=x && x<='$a') ? 0 : NaN' >> gnu-script
            echo 'set arrow from '$a',0 to '$a','$max' nohead linecolor "red"' >> gnu-script
            ;;
        * )
            echo "Error"
            ;;
    esac
    
    echo 'plot "'$1'.tmp" title "'$2'" w l linecolor "blue", f(x) title "'$3'" w l linecolor "red"' >> gnu-script

    gnuplot gnu-script;
}

TracGraphGnuplot "Tableau" "fonction d'onde" "V(x)" "" "" $borneInf $borneSup $8 $9x					

case $TypePot in
        1 )
            if ! [[ -e Nul ]]
            then
                mkdir Nul
            fi
            mv FctOnde.jpeg Nul
            mv TableauDeValeurs.txt Nul
            ;;
        2 )
            if ! [[ -e Rectangle ]]
            then
                mkdir Rectangle
            fi
            mv FctOnde.jpeg Rectangle
            mv TableauDeValeurs.txt Rectangle
            ;;
        3 )
            if ! [[ -e Marche ]]
            then
                mkdir Marche
            fi
            mv FctOnde.jpeg Marche
            mv TableauDeValeurs.txt Marche
            ;;
        4 )
            if ! [[ -e Parabole ]]
            then
                mkdir Parabole
            fi
            mv FctOnde.jpeg Parabole
            mv TableauDeValeurs.txt Parabole
            ;;
        * )
            echo "Error"
            ;;
    esac
rm Tableau.tmp
rm gnu-script
