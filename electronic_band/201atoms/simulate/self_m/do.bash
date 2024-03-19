cp ../../self/inputmm_1 all_input
cp ../../self/outputmm_1 all_output
cp ../../self/indexmm_1 all_index
cp ../../self/diagonalmm_1 all_diagonal

python 1.py 
paste pred_target0.dat all_diagonal > temp
python add.py 
