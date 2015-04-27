max=4000
for ((i=1000; i<=$max; i+=500 )) ; 
do
    python main.py $i
done
