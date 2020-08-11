executable				= #script
output                	= #logDir/#jobName.out
error                 	= #logDir/#jobName.err
log                   	= #logDir/#jobName.log
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
request_memory			= 2000M
request_disk			= 1000M
Requirements 			= (Machine != "scorpion6.spa.umn.edu") && (Machine != "spa-sl7test.spa.umn.edu") && (Machine != "zebra01.spa.umn.edu") && (Machine != "zebra02.spa.umn.edu")  && (Machine != "zebra03.spa.umn.edu")  && (Machine != "zebra04.spa.umn.edu") && (Machine != "scorpion3.spa.umn.edu") && (Machine != "scorpion4.spa.umn.edu") && (Machine != "scorpion5.spa.umn.edu")
+JobFlavour 			= "#jobflavour"
queue
