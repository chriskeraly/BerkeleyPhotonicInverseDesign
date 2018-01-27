function result = runLumericalScript(lumerical,matName,lsfName)

done=0; lumError='';
save(matName,'done','lumError','-append','-v7.3');

[status,result] = system([lumerical ' -run ' lsfName ' &']);

while(~done)
    try
        load(matName,'done','lumError');
    catch
        pause(5);
    end
    pause(1);
end

if(~isempty(lumError))
    error(lumError);
end

end
