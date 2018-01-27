
function[data] = permuteLumData(data,p,dim)
    if( length(p) == 4 )
        if( length(dim)==4 ) % many frequencies
            data = squeeze( mat2cell(permute(data,p), dim(2), dim(1), dim(3), ones(1,dim(4))) );
        else % only one frequency
            data = {permute(data,p(1:3))};
        end
    elseif( length(p) == 3)
        if( length(dim)==3 ) % many frequencies
            if(dim(1)==1 || dim(2)==1)
                data = {permute(data,p(1:3))};
            else
                data = squeeze( mat2cell(permute(data,p), dim(2), dim(1), ones(1,dim(3))) );
            end
        else % only one frequency
            data = {permute(data,p(1:2))};
        end
    end
end