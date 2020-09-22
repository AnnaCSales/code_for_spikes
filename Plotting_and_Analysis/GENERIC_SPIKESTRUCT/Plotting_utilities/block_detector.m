function [blocktimes] = block_detector(ts)
%Detects distinct blocks of data within continuous timestamps - for OEP
%recs in multiple blocks.
blocktimes=[];
if sum(diff(ts)>0.1)
    jump_inds=find(diff(ts)>0.1);
    for j=1:length(jump_inds)
        if j==1
            blocktimes(j,1)=1;
            blocktimes(j,2)=ts(jump_inds(j));
        else
            blocktimes(j,1)=ts(jump_inds(j-1)+1);
            blocktimes(j,2)=ts(jump_inds(j));
        end
    end

end

if length(jump_inds)
     fprintf('Multiple blocks detected \n')
     for j=1:length(jump_inds)
         fprintf('Block %d, %.2d to %.2d s \n', j, blocktimes(j,1),blocktimes(j,2))
     
     end
end
end

