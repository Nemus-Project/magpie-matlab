d = dir('*_help.m');
files = {d.name};  % cell array of filenames

for i = 1:numel(files)    
    publish(files{i}, 'format', 'html','outputDir', './html','evalCode',false);
end
