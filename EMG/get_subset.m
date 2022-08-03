
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 3 Mar 2021. Last revision: 17 Mar 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [structout] = get_subset (mystruct, labeltoDel)
% This function will remove the whole line(s) of a struct which carries certain label
%    mystruct: input struct
%    labeltoDel: label/name to be removed

temp = struct2cell(mystruct);

for i=size(temp,2):-1:1
    todel = strcmp(temp(1,i),labeltoDel);   % Check if the string matches to delete
    if(~todel)
        temp(:,i) = [];        % If yes then remove the whole row
    end
end

%temp
structout = cell2struct(temp,fieldnames(mystruct));  % Convert back to the struct data type
