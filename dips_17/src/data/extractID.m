function extractID(input_file, output_file)

	q_list 	= importdata(input_file);
	q_size	= length(q_list.data);
	id_list = cell(1, q_size);

	for i = 1:q_size,
		s = strsplit(q_list.textdata{i}, {'_'});
		id_list{i} = s{end}(3:end);
    end
    
    fout = fopen(output_file,'w');
    for i = 1:length(id_list),
        fprintf(fout,'%s\n',id_list{i});
    end
    fclose(fout);
    
end