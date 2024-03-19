python train.py \
	--d_model 128 \
	--h_dim 256 \
	--num_layers 3 \
	--conv_type 'cgcnn' \
	-bs 8 \
	-dr 1 \
	-e 200 \
	--optim 'adamw' \
	-lr 1e-3 \
	-wd 0 \
	--gpu_idx 0 \
	--deterministic True \
	--log_wandb \
	
