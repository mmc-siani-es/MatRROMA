function [ref_vect,pairs,r]=arrangesets(Sets,setref)
%%Function for locating reference sensors and sorting them according to reference set
%%%%Location of the reference sensors
r = max(cellfun(@max, Sets));

pairs=[];
ref_vect=[];
keepnorm=setref;
[Sets,pairs,ref_vect,keepnorm]=arrange(Sets,pairs,ref_vect,keepnorm);

while ~isempty(keepnorm)
    [Sets,pairs,ref_vect,keepnorm]=arrange(Sets,pairs,ref_vect,keepnorm);
end

if size(pairs,1)~=length(Sets)-1
    error("There is no common reference sensors between one or more Sets")
end
end

function [Sets,pairs,ref_vect,keepnorm]=arrange(Sets,pairs,ref_vect,keepnorm)
keeprow=Sets{keepnorm(1)};

Sets{keepnorm(1)}=[];
un=[];
ref_vect_temp=[];

for i=1:length(Sets)
    if isempty(Sets{i})
        continue; %Ignore visited sets
    end
     %Find the intersection of sensors
    common_sens=intersect(keeprow, Sets{i})';
    common_sens=common_sens(:)';

    if ~isempty(common_sens)
        %Avoid adding sets that have already been linked to another
        % reference
        if isempty(pairs)||~any(pairs(:,2)==i)
            un=[un;i];
            ref_vect_temp=sort(ref_vect_temp,2);
            ref_vect_temp=[ref_vect_temp; common_sens];
        end
        if (~isempty(pairs) && any(pairs(:, 2)==i)&&~all(all(ref_vect == ref_vect(1,:),2)))
            error("In roving references case a set can not be linked to more than 1 set")
        end
    end
end

%If new connected sets were found
if ~isempty(un)
    %Verify that the number of reference sensors is consistent
    if ~isempty(ref_vect) && size(ref_vect_temp, 2) ~= size(ref_vect, 2)
        error('There are different numbers of reference sensors between sets')
    end

    %Save the reference sensors
    ref_vect=[ref_vect;ref_vect_temp];

    %Register the pairs [Parent_Set, Child_Set]
    pairs_temp=zeros(length(un),2);
    pairs_temp(:,1)=keepnorm(1);
    pairs_temp(:,2)=un;

    pairs=[pairs;pairs_temp];

    %Add the new children to the queue for future iterations
    keepnorm=[keepnorm;un];
end

%Remove the current node from the work queue
keepnorm(1)=[];
end