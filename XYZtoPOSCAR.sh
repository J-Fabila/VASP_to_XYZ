if [ -f $2 ]
then
mv $2 $2other
fi
N=$(head -1 $1)
head -$(($N+2)) $1 | tail -$N |  awk '{print $1 }' | sort -u >> aux
M=$(wc -l aux | awk '{ print $1 }')
for ((i=1;i<$(($M+1));i++))
do
echo -n "$(head -$i aux | tail -1  )  " >>$2
done
echo " " >>$2
for ((i=1;i<$(($M+1));i++))
do
echo -n "$( grep "$(head -$i aux | tail -1 )" $1 | wc -l )  ">> $2
done
echo "  " >> $2
Sel=$(head -$(($N+2)) $1 | tail -$N |  awk '{print $5 }' | grep . | wc -l )
if [ $Sel -gt 0 ]
then
echo "Selective dynamics " >> $2
fi

echo "Cartesian" >> $2
rm aux
head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print $2 "  "$3 "  " $4 }' >> aux1
head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print  $5 "  " $6 "  " $7}' >> aux2
if [ $Sel -gt 0 ]
then
Nl=$(cat aux2 | wc -l )
for((i=1;i<$(($N+1)); i++))
do
cont=$(head -$i aux2 | tail -1 | tr '0' 'F' | tr '1' 'T' )
ki=$(echo $cont | wc -c)
if [ $ki -gt 1 ]
then
echo "$cont" >>aux3
else
echo " T   T   T " >> aux3
fi
done
paste aux1 aux3 >> $2
rm aux1 aux2   aux3
else 
paste aux1 aux2 >>$2
rm aux2 aux1
fi
