function [table] = InputExcel(file_path ,day)

path = file_path + day;
table = {};


list = dir(path);
files = [];


for k = 1:length(list)
    [~, name,ext] = fileparts(list(k).name);
    
    name_str = convertCharsToStrings(name);
    ext_str = convertCharsToStrings(ext);

    if strcmp( ext_str, '.xlsx')
        file_str = name_str +ext_str;

        files = [files; file_str] ;
    end

end


new_path = path + "\" + files(1);

sheets = sheetnames(new_path);

for sheet = 1:length(sheets)
    current_table = readtable(new_path, 'Sheet', sheets(sheet), 'Format', 'auto');
    table{sheet} = current_table;
end

end