clear
restoredefaultpath
addpath(genpath('submodules'));
addpath('~/git/gisogrimm/tascar/scripts')
check_for_tascar()

plugin = 'hoa2d_enc';
crop = 512;
filename='response_temp.mat'

out_dir = '~/Documents/data';
check_output_dir_exists(out_dir)


shFolder = "/home/fjm518/Desktop/20220526_215439_Thiemann2019_horiz_fir_matrix/ambix_fir_mat_Thiemann2019_horiz_o";
shSuffix = "a/ambix2array_fir.wav";

load Thiemann2019.mat
nmic = 4;
eval_az_deg = az_deg;
eval_inc_deg = inc_deg;






N = length(eval_az_deg);

orders = [3 7 23];
N = length(eval_az_deg);
response = zeros(length(orders),crop,nmic,length(eval_az_deg));
for j = 1:length(orders)
        
        order = orders(j)
        components = 2*order+1;


        %load the Thiemann SH components
        [Thie_SH,Fs] = audioread(shFolder+string(order)+shSuffix);



        x = zeros(crop,components,N); % preallocated space for rendered source irs
        for iaz=1:N

            tsc_file = sprintf('%03d.tsc',iaz);
            wav_file = sprintf('%03d.wav',iaz);
            
            %% create tascar scn file
            n_doc = tascar_xml_doc_new();
            % session (root element):
            n_session = tascar_xml_add_element(n_doc,n_doc,'session');
            n_scene = tascar_xml_add_element(n_doc,n_session,'scene');
            % source
            n_src = tascar_xml_add_element(n_doc, n_scene, 'source');
            n_sound = tascar_xml_add_element(n_doc, n_src, 'sound',[],'x','1');
            % receiver:
            n_rec = tascar_xml_add_element(n_doc,n_scene,'receiver',[],...
                'type', plugin,...
                'order',int2str(order),...
                'name', sprintf('out_%d',iaz));
            tascar_xml_add_element(n_doc, n_rec, ...
                                   'orientation',sprintf('0 %g 0 0',eval_az_deg(iaz)));

            tascar_xml_save( n_doc, tsc_file);

            system(['LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH ',...
                    sprintf('tascar_renderir %s --srate %d -o %s', ...
                            tsc_file, Fs, wav_file)]);
            
            [xtmp,~] = v_readwav(wav_file);
            x(:,:,iaz) = xtmp(1:crop,:,:);
            
            delete(tsc_file,wav_file);
        end 
        


        y = [];
        for imic = nmic:-1:1
            SH_yfull = reshape(squeeze(Thie_SH(:,imic)),[],(order+1)^2);
            SH_y =SH_yfull(:,any(SH_yfull));
            for idoa = N:-1:1
                y(:,imic,idoa) = sum(fftfilt(SH_y,x(:,:,idoa)),2);
            end
        end

        response(j,:,:,:) = y;
end
spkno = 2.*(orders)+1;
save(filename,'response','orders','spkno')
