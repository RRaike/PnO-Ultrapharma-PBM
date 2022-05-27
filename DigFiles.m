function [result, names] = DigFiles(path)

list = FindAllInFolder(path);
result = [];
names = [];

if ~isempty(list)
    for k = 1: length(list)
    list_element = [list{k}];
    list_string = convertCharsToStrings(list_element(:));
    new_path = path + '\'+ list_string;

    result = [result; new_path];
    names = [names; list_string];
    end
end
end