function [list] = FindAllInFolder(path)

list = {};

list_1 = dir(path);

for k = 1:length(list_1)
    
    if ~strcmp(list_1(k).name, '.') && ~strcmp(list_1(k).name, '..')
        if list_1(k).bytes == 0
            
            list{end+1} = list_1(k).name;
       
        end
    end
end

end