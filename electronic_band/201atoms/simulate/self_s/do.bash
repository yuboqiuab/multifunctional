cp ../../self/inputss_1 all_input
cp ../../self/outputss_1 all_output
cp ../../self/indexss_1 all_index
cp ../../self/diagonalss_1 all_diagonal
python 1.py 
paste pred_target0.dat all_diagonal > see
python add.py 
