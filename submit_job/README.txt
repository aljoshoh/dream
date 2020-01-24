# Usage #

In order to get some info about usage, run "execute_calc_slurm_multi.sh" !
As example:

./tools/submit_job/execute_calc_slurm.sh 1 3 test_script.R tools/submit_job rstudio_server_icb 2

                        		 ^
					 index of first parallel job
					   ^
					   index of last parallel job
			       		     ^
			       		     name of script should be executed
								 ^
								 directory where the script lays
									    ^
									    docker image in the image folder (image folder is hardcoded)
											       ^
											       tasks per array job
