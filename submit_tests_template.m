addpath INSERT_BASE_DIR_HERE 
cd INSERT_INSTANCE_DIR_HERE
[rm_old,relgap,rm_new,cuts,rv] = liftedRLT_fun('INSERT_FILE_NAME_HERE');
save('INSERT_RESULTS_DIR_HERE/INSERT_FILE_NAME_HERE');