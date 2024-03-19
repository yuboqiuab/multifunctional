rm all_force.dat

cd step0
cd sinput/
python t1.py
cp input_1 output_1 ../sforce/
cd ../moinput/
python t1.py
cp output_1 input_1 ../moforce/
cd ../sforce/
python 1.py --resume
cd ../moforce/
python 1.py --resume
cd ../
cd relax/
bash do.bash
cd ../
cat relax/maxf >> ../all_force.dat
cd ../

for ((j=1;j<=275;j=j+1))
    do
    i=$(echo "("$j"-1)" | bc -l)
    cp -r step"$i" step"$j"
    cd step"$j"
    bash all.bash
    cat relax/maxf >> ../all_force.dat  
    cd ../
    done
