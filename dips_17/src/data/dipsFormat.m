function data = dipsFormat(opt)
% To transform separate correlation files into DIPS input format

	% define dir of correlation matrices
	if isfield(opt, 'corr_dir')
		corr_dir = opt.corr_dir;
	else
		err('... please provide the directory of correlation matrices ...');
	end
	fstruct = dir([corr_dir, '/*.mat']);
	snap_num = length(fstruct);

	% define #node and #edge
	load([corr_dir, '/', fstruct(1).name]);
	node_num = size(frame.mat,1);
	edge_num = node_num*(node_num-1)/2;

	% create DIPS input structure 
	data = [];
	data.X 		= zeros(snap_num, edge_num);
	data.label	= zeros(snap_num, 1);
	data.site 	= cell(1, snap_num);
	data.id		= cell(1, snap_num);
	if isfield(opt, 'WL')
		data.init	= zeros(snap_num, 1);
		data.end 	= zeros(snap_num, 1);
	end

	% transform matrices to feature vectors
	for i = 1:snap_num
		% load time-series of each subject
        disp(i);
	    load([corr_dir, '/', fstruct(i).name]);
	    data.X(i,:) 	= mat2vec(frame.mat)';
	    data.site{i}	= frame.site;
		data.id{i}		= frame.id;
		if isfield(opt, 'WL')
			data.init(i,1)	= frame.init;
			data.end(i,1) 	= frame.end;
		end
	end

	% add label to each functional correlation matrix
	opt.meta_dir = [corr_dir,'/../pheno'];
	mkdir_if_not_exist(opt.meta_dir)
	data.gnd = loadPheno(data.id, 8, opt)';

	% add network adjacent matrix
	data.W = sparse(genLineGraph(node_num));

	% save dataset to mat file
	save([corr_dir,'/../dips_data.mat'],'data');

end