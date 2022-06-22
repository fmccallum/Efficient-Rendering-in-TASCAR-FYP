% script to show use of this tool

%% make sure we aren't using anything we shouldn't be
clear
restoredefaultpath
addpath(genpath('submodules'));
addpath('~/git/gisogrimm/tascar/scripts')
check_for_tascar()


plugin = 'vbap_fredlut';
filename='response_vbaplut.mat'

out_dir = '~/Documents/data';
check_output_dir_exists(out_dir)


spk_file = 'temp.spk';


load Thiemann2019.mat

eval_az_deg = az_deg;
spkno = [16,30,64,120]; %number of virtual speakers


iarray =  @Sampled_20220225_Thiemann2019_bte_horiz;

%%
vspk = struct();
vspk.ema_fcn = iarray;
vspk.fs = fs;


%Load in ground truth data

response = zeros(length(spkno),512,4,length(eval_az_deg));
for j = 1:length(increments)
        spkno(j)
        %Set virtual speakers
        az_deg = (-180:increments(j):179).';
        inc_deg = 90*ones(size(az_deg));
        r = ones(size(az_deg));
        vspk.doa_unit_vec = mysph2cart(pi/180*az_deg,...
                                   pi/180*inc_deg,...
                                   r);
        [vspk_ir, cart_vec] = fcn_20220203_01_get_vspk_to_array_irs(vspk);
        [nsamples,nmic,nvspk] = size(vspk_ir);


    % convert vspk positions to TASCAR's spherical coordinates
            tascarSphCoord.r_m = sqrt(sum(cart_vec.^2,2));
            tascarSphCoord.el_deg = 90 - (180/pi) * acos(cart_vec(:,3)./tascarSphCoord.r_m);
            tascarSphCoord.az_deg = (180/pi) * atan2(cart_vec(:,2),cart_vec(:,1));

    % write an spk file which defines (virtual) speaker positions
%     fcn_20200129_02_write_spk_file(fullfile(temp_dir,spk_file), tascarSphCoord);
            fcn_20200129_02_write_spk_file(spk_file, tascarSphCoord);



        N = length(eval_az_deg);
        x = zeros(512,nvspk,N); % preallocated space for rendered source irs
        for iaz=1:N

            tsc_file = sprintf('%03d.tsc',iaz);
            wav_file = sprintf('%03d.wav',iaz);
            
            %% create tascar scn file
            n_doc = tascar_xml_doc_new();

            n_session = tascar_xml_add_element(n_doc,n_doc,'session');
            n_scene = tascar_xml_add_element(n_doc,n_session,'scene');
            % source
            n_src = tascar_xml_add_element(n_doc, n_scene, 'source');
            n_sound = tascar_xml_add_element(n_doc, n_src, 'sound',[],'x','1');
            % receiver:
            n_rec = tascar_xml_add_element(n_doc,n_scene,'receiver',[],...
                         'type', plugin,...
                         'name', sprintf('out_%d',iaz),...
                         'layout', spk_file);
            tascar_xml_add_element(n_doc, n_rec, ...
                                   'orientation',sprintf('0 %g 0 0',eval_az_deg(iaz)));

            tascar_xml_save( n_doc, tsc_file);

            %% execute tascar to obtain rendered impulse responses
            system(['LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ',...
                    sprintf('tascar_renderir %s --srate %d -o %s', ...
                            tsc_file, vspk.fs, wav_file)]);
            %% read in the result
            [xtmp,~] = v_readwav(wav_file);
            x(:,:,iaz) = xtmp(1:512,:,:);
            
            delete(tsc_file,wav_file);
        end 

        xcrop = x(1:512,:,:);

        %% compute the per mic signals
        y = [];
        for idoa = N:-1:1
            for imic = nmic:-1:1
                y(:,imic,idoa) = sum(fftfilt(permute(vspk_ir(:,imic,:),[1 3 2]),...
                                             xcrop(:,:,idoa)),...
                                     2);
            end


        end

        response(j,:,:,:) = y;
end

save(filename,'response','spkno')
