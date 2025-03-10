# Ensure we don't have any conflicting modules loaded
ml purge; 

#Insert container
container=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_gaswar/gaswar_container.sif;

# You can launch jupyter notebook or lab, but you must specify the config file as below: 
apptainer exec $container jupyter lab --config="${CONFIG_FILE}";
