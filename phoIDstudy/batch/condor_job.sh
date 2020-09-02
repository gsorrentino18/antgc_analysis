executable				= #script
output                	= #logDir/#jobName.out
error                 	= #logDir/#jobName.err
log                   	= #logDir/#jobName.log
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
request_memory			= 2000M
request_disk			= 1000M
#Requirements 			=  (Machine != "scorpion10.spa.umn.edu") && (Machine != "scorpion6.spa.umn.edu") &&(Machine != "scorpion1.spa.umn.edu") && (Machine != "spa-sl7test.spa.umn.edu") && (Machine != "zebra01.spa.umn.edu") && (Machine != "zebra02.spa.umn.edu")  && (Machine != "zebra03.spa.umn.edu")  && (Machine != "zebra04.spa.umn.edu") && (Machine != "scorpion3.spa.umn.edu") && (Machine != "scorpion4.spa.umn.edu") && (Machine != "scorpion5.spa.umn.edu")
# Requirements 			=  (Machine == "scorpion3.spa.umn.edu") || (Machine == "scorpion4.spa.umn.edu") ||(Machine == "scorpion5.spa.umn.edu") ||(Machine == "scorpion5.spa.umn.edu") ||(Machine == "scorpion7.spa.umn.edu") ||(Machine == "scorpion8.spa.umn.edu") ||(Machine == "scorpion9.spa.umn.edu") ||(Machine == "scorpion10.spa.umn.edu")
#Requirements 			=  (Machine == "scorpion6.spa.umn.edu")
Requirements 			=  (Machine == "scorpion3.spa.umn.edu") || (Machine == "scorpion4.spa.umn.edu") || (Machine == "scorpion5.spa.umn.edu") || (Machine == "scorpion6.spa.umn.edu") || (Machine == "scorpion7.spa.umn.edu") || (Machine == "scorpion8.spa.umn.edu") || (Machine == "scorpion9.spa.umn.edu") || (Machine == "scorpion10.spa.umn.edu") || (Machine == "scorpion11.spa.umn.edu") || (Machine == "scorpion12.spa.umn.edu") || (Machine == "scorpion13.spa.umn.edu") || (Machine == "scorpion14.spa.umn.edu") || (Machine == "scorpion15.spa.umn.edu") || (Machine == "scorpion16.spa.umn.edu") || (Machine == "scorpion17.spa.umn.edu") || (Machine == "scorpion18.spa.umn.edu") || (Machine == "scorpion19.spa.umn.edu") || (Machine == "scorpion20.spa.umn.edu") || (Machine == "scorpion21.spa.umn.edu") || (Machine == "scorpion22.spa.umn.edu") || (Machine == "scorpion23.spa.umn.edu") || (Machine == "scorpion24.spa.umn.edu") || (Machine == "scorpion25.spa.umn.edu") || (Machine == "scorpion26.spa.umn.edu") || (Machine == "scorpion27.spa.umn.edu") || (Machine == "scorpion28.spa.umn.edu") || (Machine == "scorpion29.spa.umn.edu") || (Machine == "scorpion30.spa.umn.edu") || (Machine == "scorpion31.spa.umn.edu") || (Machine == "scorpion32.spa.umn.edu") || (Machine == "scorpion33.spa.umn.edu") || (Machine == "scorpion34.spa.umn.edu") || (Machine == "scorpion35.spa.umn.edu") || (Machine == "scorpion36.spa.umn.edu")
#next_job_start_delay	= 3
+JobFlavour 			= "#jobflavour"
queue
