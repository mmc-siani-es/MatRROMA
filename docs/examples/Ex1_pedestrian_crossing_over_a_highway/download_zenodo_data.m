function download_zenodo_data()

% Gets the folder where the script is located
script_dir=fileparts(mfilename('fullpath'));
url_zenodo='https://zenodo.org/records/21703314/files/Experimental_data_Elevated_pedastrian_crossing.zip?download=1';
zip_filename=fullfile(script_dir,'Experimental_data_Elevated_pedastrian_crossing');
required_files = {'Lvl_crossing_P1_P2.dat','Lvl_crossing_P2_P3.dat','Lvl_crossing_P3_P4.dat','Lvl_crossing_P4_P5.dat','Lvl_crossing_P5_P6.dat'};

% Check if ALL 6 required files exist
needs_download = false;
for i = 1:length(required_files)
    if ~isfile(fullfile(script_dir, required_files{i}))
        needs_download = true;
        break;
    end
end

if needs_download
    disp('Downloading required data from Zenodo...');
    try
        websave(zip_filename, url_zenodo);

        % 1. Create a safe temporary directory
        safe_temp_dir=fullfile(tempdir, 'zenodo_extraction_zone');
        if isfolder(safe_temp_dir)
            rmdir(safe_temp_dir, 's');
        end
        mkdir(safe_temp_dir);


        disp('Extracting files safely...');
        extracted_paths=unzip(zip_filename, safe_temp_dir);

        % 3. Move ONLY the 6 required files to the project folder
        disp('Filtering and moving only the necessary files...');
        for k = 1:length(extracted_paths)
            if isfile(extracted_paths{k})
                [~,name,ext]=fileparts(extracted_paths{k});
                file_only=[name,ext];

                % Check if this specific file is in your list of 6
                if any(strcmp(required_files, file_only))
                    dest_path=fullfile(script_dir, file_only);
                    movefile(extracted_paths{k},dest_path,'f');
                end
            end
        end

        % 4. Clean up the temporary zone 
        rmdir(safe_temp_dir,'s');
        if isfile(zip_filename)
            delete(zip_filename);
        end
        disp('Data downloaded and extracted successfully!');

    catch ME
        disp(['Error during download/extraction: ',ME.message]);
    end
end
end