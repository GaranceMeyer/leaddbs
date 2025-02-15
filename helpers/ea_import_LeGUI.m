function ea_import_localization(dataset, subjId, LeGUIdir, space)
% Import localization from other tool into Lead-DBS subject folder.
arguments
    dataset  {mustBeTextScalar} % Dataset path
    subjId   {mustBeTextScalar} % Subject ID
    LeGUIdir {mustBeTextScalar} % Path of LeGUI subject directory
    space    {mustBeMember(space, {'native', 'mni'})} = 'mni' % Space of the coordinates

    % Load LeGUI coordinates
    load([LeGUIdir filesep 'Registered' filesep 'Electrodes.mat']);

    % Identify electrodes based on contact labels
    c_labels = cellfun(@(x) regexprep(x, '\d+$', ''), ElecMapRaw(:,1), 'UniformOutput', false);
    c_numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), ElecMapRaw(:,1));
    el_names = unique(c_labels(DepthElecRaw)); % here, only select electrodes that are depth electrodes
    fprintf('Processing subject %s - Found %d electrodes\n', subjID, length(el_names))

    % Get skull voxel coordinates for detection of most distal contacts
    V=spm_vol([LeGUIdir filesep 'Registered' filesep 'MRBone.nii']);
    [Y,XYZ] = spm_read_vols(V); % read the image and mmcoordinates
    skull_coords = XYZ(:, (Y(:)>0.5));

    for ii = 1:length(el_names)
        fprintf('Processing electrode %s.\n', el_names{ii})

        reco.props(ii).elmodel = 'Multiple';
        reco.props(ii).manually_corrected = 1;

        % find indices of contacts corresponding to this electrode
        el_idx = find(strcmp(c_labels, el_names{ii}));

        % reorder contacts 
        native_coords = ElecXYZRaw(el_idx, :);
        proj_coords   = ElecXYZProjRaw(el_idx, :);
        mni_coords    = ElecXYZMNIProjRaw(el_idx, :);
        
        % will arbitrarily start from largest |X|, or if electrodes are more vertical,
        % from largest z and then will order contacts based on distance 
        if abs(max(native_coords(:,1))-min(native_coords(:,1)))>10
            [~, start_idx] = max(abs(native_coords(:,1)));
        else 
            [~, start_idx] = max(native_coords(:,3));
        end

        d = sqrt(sum((native_coords - native_coords(start_idx, :)).^2, 2));
        [~,sort_idx] = sort(d, 'descend'); 

        % then check whether the first or last contact is closest to
        % skull and reorder if contacts are ordered in the wrong direction
        if min(sqrt(sum((skull_coords' - native_coords(sort_idx(1),:)).^2,2))) < ... 
                min(sqrt(sum((skull_coords' - native_coords(sort_idx(end),:)).^2,2))) 
            [~,sort_idx] = sort(d, 'ascend'); 
        end
        
        % fill coordinates in reco
        reco.native.coords_mm{ii} = native_coords(sort_idx, :); 
        reco.scrf.coords_mm{ii}   = proj_coords(sort_idx, :); 
        reco.mni.coords_mm{ii}    = mni_coords(sort_idx, :); 

        if ~issorted(c_numbers(el_idx(sort_idx)),'ascend')
            ea_warning('Labels are inconsistent with automatic contact ordering.')
        end 
        
        % check - pinker contacts should more distal
        for ci = 1:length(reco.mni.coords_mm{ii})
            ea_plotsphere(reco.mni.coords_mm{ii}(ci, :), 2, [0.05*ci, 0, 1-(0.05*ci)]); 
        end 

        % electrode model
        
        models = ea_resolve_elspec;
        avg_spacing = mean(sqrt(sum(diff(native_coords(sort_idx, :)).^2, 2)));
        prompt = {'Pick electrode model - ', sprintf('Mean spacing: %.2f mm,',avg_spacing), sprintf('%d contacts',size(native_coords,1)) };
        index = listdlg('PromptString', prompt, 'ListString', models, 'SelectionMode', 'single', 'CancelString', 'Cancel');
        if ~isempty(index)
            elmodel = models{index};
        else
            return;
        end

        options.elmodel = elmodel;
        options = ea_resolve_elspec(options);
        
        reco.props(ii).multiple_elmodel = elmodel; 
        reco.props(ii).elname = el_names{ii}; 
        reco.props(ii).labels = ElecMapRaw((el_idx(sort_idx)),1); 
        
        % fill head, tail and trajectory - note that this is far suboptimal
        % for SEEG, to be fixed later 
        reco.native.markers(ii).head = reco.native.coords_mm{ii}(1, :);
        reco.native.markers(ii).tail = reco.native.coords_mm{ii}(4, :);  

        [xunitv, yunitv] = ea_calcxy(reco.native.markers(ii).head, reco.native.markers(ii).tail);
        reco.native.markers(ii).x = reco.native.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
        reco.native.markers(ii).y = reco.native.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
        [~, reco.native.trajectory(ii)] = ea_resolvecoords(reco.native.markers(ii), elmodel);

        reco.scrf.markers(ii).head = reco.scrf.coords_mm{ii}(1, :);
        reco.scrf.markers(ii).tail = reco.scrf.coords_mm{ii}(4, :);  

        [xunitv, yunitv] = ea_calcxy(reco.scrf.markers(ii).head, reco.scrf.markers(ii).tail);
        reco.scrf.markers(ii).x = reco.scrf.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
        reco.scrf.markers(ii).y = reco.scrf.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
        [~, reco.scrf.trajectory(ii)] = ea_resolvecoords(reco.scrf.markers(ii), elmodel);

        reco.mni.markers(ii).head = reco.mni.coords_mm{ii}(1, :);
        reco.mni.markers(ii).tail = reco.mni.coords_mm{ii}(4, :);  

        [xunitv, yunitv] = ea_calcxy(reco.mni.markers(ii).head, reco.mni.markers(ii).tail);
        reco.mni.markers(ii).x = reco.mni.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
        reco.mni.markers(ii).y = reco.mni.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
        [~, reco.mni.trajectory(ii)] = ea_resolvecoords(reco.mni.markers(ii), elmodel); 

    end % for ii = 1:length(el_names)

    % save reco
    derivativesFolder = fullfile(dataset, 'derivatives', 'leaddbs', ['sub-', subjId]);
    rawFolder = fullfile(dataset, 'rawdata', ['sub-', subjId]);
    ea_mkdir({derivativesFolder; rawFolder});
    bids = BIDSFetcher(dataset);
    subj = bids.getSubj(erase(subjID, textBoundary("start") + 'sub-'));

    ea_mkdir(fileparts(subj.recon.recon));
    save(subj.recon.recon, 'reco');

    % import other files
    % CT, MRI
    ea_mkdir(subj.coregDir);
    copyfile(fullfile(LeGUIdir, 'registered', 'CT.nii'), ...
        );
    copyfile(fullfile(LeGUIdir, 'registered', 'MR.nii'), ...
        );

    % segmentations
    copyfile(fullfile(LeGUIdir, 'registered', 'MRCSF.nii'), ...
        );
    copyfile(fullfile(LeGUIdir, 'registered', 'MRGray.nii'), ...
    );
    copyfile(fullfile(LeGUIdir, 'registered', 'MRWhite.nii'), ...
        );

    % normalization matrices
    ea_mkdir(subj.normDir);

    copyfile(fullfile(LeGUIdir, 'registered', 'iy_MR.nii'), ...
        );
    copyfile(fullfile(LeGUIdir, 'registered', 'y_MR.nii'), ...
        );

    % brain surface
    ea_mkdir(fullfile(subj.subjDir, 'surface'));
    copyfile(fullfile(LeGUIdir, 'registered', 'Surfaces.mat'), ...
        );

    % copy all files but these in miscellaneous


end
