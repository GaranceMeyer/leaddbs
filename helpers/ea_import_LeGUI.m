function ea_import_localization(dataset, subjID, LeGUIdir, space)
% Import localization from other tool into Lead-DBS subject folder.
arguments
    dataset  {mustBeTextScalar} % Dataset path
    subjID   {mustBeTextScalar} % Subject ID
    LeGUIdir {mustBeTextScalar} % Path of LeGUI subject directory
    space    {mustBeMember(space, {'native', 'mni'})} = 'mni' % Space of the coordinates

    % Load LeGUI coordinates
    load([LeGUIdir filesep 'Registered' filesep 'Electrodes.mat']);

    % Identify electrodes based on contact labels
    c_labels = cellfun(@(x) regexprep(x, '\d+$', ''), ElecMapRaw(:,1), 'UniformOutput', false);
    c_numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), ElecMapRaw(:,1));
    el_names = unique(c_labels);
    fprintf('Processing subject %s - Found %d electrodes\n', subjID, length(el_names))

    % Get skull voxel coordinates for detection of most distal contacts
    V=spm_vol([LeGUIdir filesep 'Registered' filesep 'MRBone.nii']);
    [Y,XYZ] = spm_read_vols(V); % read the image and mmcoordinates
    skull_coords = XYZ(:, (Y(:)>0.5));

    for ii = 1:length(el_names)
        fprintf('Processing electrode %s.\n', el_names{ii})

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
        [~,sort_idx] = sort(d, 'ascend'); 

        % then check whether the first or last contact is closest to
        % skull and reorder if contacts are ordered in the wrong direction
        if min(sqrt(sum((skull_coords' - native_coords(sort_idx(1),:)).^2,2))) > ... 
                min(sqrt(sum((skull_coords' - native_coords(sort_idx(end),:)).^2,2))) 
            [~,sort_idx] = sort(d, 'descend'); 
        end
        
        % fill coordinates in reco
        reco.native.coords_mm{ii} = native_coords(sort_idx, :); 
        reco.scrf.coords_mm{ii}   = proj_coords(sort_idx, :); 
        reco.mni.coords_mm{ii}    = mni_coords(sort_idx, :); 

        if ~issorted(c_numbers(el_idx(sort_idx)),'descend')
            ea_warning('Labels are inconsistent with automatic contact ordering.')
        end 

%         for ci = 1:length(reco.mni.coords_mm{ii})
%             ea_plotsphere(reco.mni.coords_mm{ii}(ci, :), 2, [0.05*ci, 0, 1-(0.05*ci)]); 
%         end 
        
        models = ea_resolve_elspec;
        avg_spacing = mean(sqrt(sum(diff(native_coords(sort_idx, :)).^2, 2)));
        prompt = {'Pick electrode model - ', sprintf('Mean spacing: %.2f mm,',avg_spacing), sprintf('%d contacts',size(native_coords,1)) };
        index = listdlg('PromptString', prompt, 'ListString', models, 'SelectionMode', 'single', 'CancelString', 'Cancel');
        if ~isempty(index)
            elmodel = models{index};
        else
            return;
        end


    end 









    % fill coords and head / tail

    % save reco
    [~, reco.(space).trajectory] = ea_resolvecoords(reco.(space).markers, elmodel);

    bids = BIDSFetcher(dataset);
    subj = bids.getSubj(erase(subjID, textBoundary("start") + 'sub-'));

    ea_mkdir(fileparts(subj.recon.recon));
    save(subj.recon.recon, 'reco');

    % import other files
    % CT, MRI

    % segmentations

    % normalization matrices

    % brain surface

    % copy all files but these in miscellaneous


end
