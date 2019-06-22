#! /bin/bash
N=$(head -1 $1)
head -$(($N+2))  $1 | tail -$N |  awk '{print $1 }' | sort -u >> aux
M=$(wc -l aux | awk '{ print $1 }')
for ((i=1;i<$(($M+1));i++))
do
echo -n "$(head -$i aux | tail -1  )  " >> output_file
done
echo " " >> output_file
for ((i=1;i<$(($M+1));i++))
do
echo -n "$( grep -w "$(head -$i aux | tail -1 )" $1 | wc -l )  ">> output_file
done
echo "  " >> output_file
Sel=$(head -$(($N+2)) $1 | tail -$N |  awk '{print $5 }' | grep . | wc -l )
if [ $Sel -gt 0 ]
then
echo "Selective dynamics " >> output_file
fi

echo "Cartesian" >> output_file
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
paste aux1 aux3 >> output_file
rm aux1 aux2   aux3
else
paste aux1 aux2 >> output_file
rm aux2 aux1
fi
if [ $# -gt 1 ]
then
   if [ -f $2 ]
   then
      mv $2 other_$2
   fi
   mv  output_file  $2
else
   cat  output_file
   rm  output_file
fi
