cp relax/new dataset/MoS2_96/structures/structures_1 
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
cd ../relax/
bash do.bash 
cd ../
