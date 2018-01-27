function done = runParallelLum(fspName, queueName)

folderPath = pwd();
simFile = [folderPath,'/',fspName];

jobFile = [folderPath,'/jobFile'];

fid = fopen(jobFile,'wt');
%fprintf(fid,['/data/lumerical/fdtd/bin/fdtd-engine-mpich2nem /data/lumerical/fdtd/bin/fdtd-engine-mpich2nem ',simFile]);
fprintf(fid,['/usr/mpi/gcc/openmpi-1.4.5/bin/mpirun -bynode /data/lumerical/fdtd/bin/fdtd-engine-ompi-lcl ',simFile]);
fclose(fid);

% Remove old log file
system(['rm -f ' strtok(fspName,'.') '_p0.log']);

% Submit job to queue
if(strcmp(queueName,'gauss'))
    system(['qsub ',jobFile,' -q ',queueName,' -N lumericalJob -l nodes=1:ppn=16,walltime=06:00:00 -u root']);
else
    system(['qsub ',jobFile,' -q ',queueName,' -N lumericalJob -l nodes=4:ppn=16,walltime=06:00:00 -u root']);
end

 tic;
 done=0;
 notdone=1;
 while( (toc<6*3600) && ~(done || ~notdone) )
     pause(1);
     [done, msg] = system(['qstat | grep ',queueName]);
     [notdone, msg] = system(['grep "Simulation completed successfully" ' strtok(fspName,'.') '_p0.log']);
 end

end
