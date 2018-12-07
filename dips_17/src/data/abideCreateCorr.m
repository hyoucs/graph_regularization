function abideCreateCorr(opt)

% = = = = = = Load window length and working directories

    if isfield(opt, 'WL')
        WL = opt.WL;    % ---- unit: TR
    end

    if isfield(opt, 'dir_ts')
        dir_ts = opt.dir_ts;
    else
        dir_ts = 'download';
    end

    if isfield(opt, 'dir_output')
        dir_output = opt.dir_output;
    else
        dir_output = 'correlation';
    end

    if isfield(opt, 'dir_prefix')
        dir_prefix = opt.dir_prefix;
    else
        if isfield(opt, 'WL')
            dir_prefix = ['wcorr-WL',num2str(WL),'-pc'];
        else
            dir_prefix = 'wcorr-session-pc';
        end
    end


% = = = = = = filt_global: cc_200, # subjects: 292
	
	% load 1-dimensional time series of selected/qualified subjects
	data = load1D(dir_ts);

	% create output directorty
	mkdir_if_not_exist(dir_output);
	mkdir_if_not_exist([dir_output,'/',dir_prefix]);
	mkdir_if_not_exist([dir_output,'/',dir_prefix,'/mat']);
	mkdir_if_not_exist([dir_output,'/',dir_prefix,'/csv']);

	% For each session, compute Pearson's correlation:
	fnum = length(data);
    for i = 1:fnum

        frame       = [];
        frame.id    = data{i}.id;
        frame.site  = data{i}.site;

        % compute correlation with sliding window
        if isfield(opt, 'WL')

            % load time-series
            [npts, nroi] = size(data{i}.data);  % npts differs from subjects
            frame_num = floor(npts/WL);    		% separate into frames of certain window length
            disp(['Subject: ', num2str(i), ' FileName:', data{i}.name]);

            % compute correlation matrices of each sliding window
            for j = 1:frame_num

                pts_beg = (j-1)*WL+1;
                pts_end = j*WL;

                frame.init = pts_beg;
                frame.end  = pts_end;
                frame.mat  = corrcoef(data{i}.data(pts_beg:pts_end,:));

                % save correlation matrices into separate mat files
                save([dir_output,'/',dir_prefix,'/mat/',num2str(frame.id),...
                	'_init_',num2str(frame.init),'.mat'], 'frame');
                dlmwrite([dir_output,'/',dir_prefix,'/csv/',num2str(frame.id),...
                	'_init_',num2str(frame.init),'.csv'], 'frame.mat');
                disp(['Created Frame: ',num2str(i),' id:',num2str(frame.id),...
                	' init:', num2str(frame.init)]);

            end

        % compute correlation with full session
        else

            frame.mat = corrcoef(data{i}.data);
            % save correlation matrices into separate mat files
            save([dir_output,'/',dir_prefix,'/mat/',num2str(frame.id),...
                '_session.mat'], 'frame');
            dlmwrite([dir_output,'/',dir_prefix,'/csv/',num2str(frame.id),...
                '_session.csv'], 'frame.mat');
            disp(['Created Frame: ',num2str(i),' id:',num2str(frame.id),...
                ' full session']);

        end 

    end

end