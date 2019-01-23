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
head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print  $5 "  " $6 "  " $7}' | tr '0' 'F' >> aux2 #tr '1' 'F'
#Si es necesario agregar una funcion que lea linea por linea, que determine si es vacia y en ese casi rellene con T
paste aux1 aux2 >> $2
rm aux1 aux2
