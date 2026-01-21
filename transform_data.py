"""
The purpose of this code is to convert the data from 2D to 3D. The data can be either charge, field or raw data.
The code is composed of a class containing the methods to convert the data. The class is called ProcessData.
The class takes as input the path to the simulation folder, the type of data to be converted, the output path,
the species, the field, the mode, the domain and the timestep. Each of these parameters will be passed by the user
and each have a default value. The default values are given in the argparser function.
the class ProcessData has three main methods one for each type of data. The methods are called convert_charge, convert_fields and
convert_raw. To run the code, the user needs to specify the path to the simulation folder and the type of data to be converted.
The user can also specify the output path, the species, the field, the mode, the domain and the timestep as flags in the command line.
an example would be:

python convert_data.py /path/to/simulation/data charge -o /path/to/output/ -s electrons -f e2 -m 1 -d 0 1 0 1 -t 100

The code will then convert the charge data from 2D to 3D and output the converted files in the specified output path.
For convenience, one can used VISX to visualize the converted files. It, however, can also be used with Python and other visualization tools.

The authos of the code are:
    - Lucas Inigo Gamiz
    - Bertrand Martinez
    - Oscar Amaro

Also, I would like to acknowledge Rafael Russo Almeida for his help.

"""



import argparse
import os
import sys
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import h5py
import shutil
import timeit
import multiprocessing as mp
from functools import partial



# Global function for parallel processing (needs to be at module level for pickling)
def process_single_dump_charge(args):
    """Process a single dump for charge conversion - designed for parallel execution"""
    dump, file_dict, dest_directory, species, i_file, mode, step_z, step_r, domain = args
    import utils
    uts = utils.utilities()
    uts.set_step_z(step_z)
    uts.set_step_r(step_r)
    if domain is not None:
        domain_dict = {}
        domain_dict['r'] = domain[:2]
        domain_dict['z'] = domain[2:]
        uts.set_domain_size(domain_dict)
    else:
        uts.set_domain_size(domain)
    
    print(f'Converting dump {i_file} (timestep {dump})')
    try:
        uts.convert_and_write_hdf5_file_densities(file_dict, dest_directory, species, i_file, mode)
        return (dump, True, None)
    except Exception as e:
        return (dump, False, str(e))


def process_single_dump_field(args):
    """Process a single dump for field conversion - designed for parallel execution"""
    dump, file_dict, dest_directory, conversion_map, mode, i_file, s_avg, step_z, step_r, domain, field_type, key = args
    import utils
    uts = utils.utilities()
    uts.set_step_z(step_z)
    uts.set_step_r(step_r)
    if domain is not None:
        domain_dict = {}
        domain_dict['r'] = [domain[0], domain[1]]
        domain_dict['z'] = [domain[2], domain[-1]]
        uts.set_domain_size(domain_dict)
    else:
        uts.set_domain_size(domain)
    
    print(f'Converting dump {i_file} (timestep {dump})')
    try:
        # Get direction value from conversion map
        dir_val = conversion_map.get(key, '3')  # Default to '3' if not found
        
        if field_type == 'e':
            uts.convert_and_write_hdf5_file_fields(file_dict, dest_directory,
                                                   f'e1_cyl_m{s_avg}', 
                                                   f'e2_cyl_m{s_avg}', 
                                                   f'e3_cyl_m{s_avg}', 
                                                   i_file, mode, dir_val)
        else:
            uts.convert_and_write_hdf5_file_fields(file_dict, dest_directory,
                                                   f'b1_cyl_m{s_avg}', 
                                                   f'b2_cyl_m{s_avg}', 
                                                   f'b3_cyl_m{s_avg}', 
                                                   i_file, mode, dir_val)
        return (dump, True, None)
    except Exception as e:
        return (dump, False, str(e))


def process_single_dump_raw(args):
    """Process a single dump for raw conversion - designed for parallel execution"""
    file_path, dest_directory, species, i_file, step_z, step_r, domain = args
    import utils
    uts = utils.utilities()
    uts.set_step_z(step_z)
    uts.set_step_r(step_r)
    if domain is not None:
        domain_dict = {}
        domain_dict['r'] = domain[:2]
        domain_dict['z'] = domain[2:]
        uts.set_domain_size(domain_dict)
    else:
        uts.set_domain_size(domain)
    
    print(f'Converting dump {i_file}')
    try:
        uts.convert_and_write_hdf5_file_raws(file_path, dest_directory, species, i_file, 
                                            if_conv_ene_to_gev=True, if_reduce=False, step=1)
        return (i_file, True, None)
    except Exception as e:
        return (i_file, False, str(e))


class ProcessData():
    def __init__(self,
                file_path : str,
                data_type: str,
                output_path : str,
                species : str,
                field : str,
                mode : int, 
                domain : list[float],
                timestep : int,
                spatial_average : bool,
                lineout : bool,
                step_z : int,
                step_r : int,
                parallel : bool = False,
                n_workers : int = None):
        
        self.file_path = file_path
        self.data_type = data_type
        self.output_path = os.getcwd() if output_path == '.' else output_path
        self.species = species
        self.field = field
        self.mode = mode
        self.domain = domain
        self.timestep = timestep
        self.spatial_average = spatial_average
        self.lineout = lineout
        self.step_z = step_z
        self.step_r = step_r
        self.parallel = parallel
        self.n_workers = n_workers if n_workers is not None else min(mp.cpu_count(), 4)
        import utils
        self.utils = utils.utilities()  # Create an instance of utilities



    def get_timestep_of_simulation(self, path_of_file) -> int:
        """
        Get the timestep of the simulation from the file name

        """
        # Use rsplit for better performance when splitting from the end
        if isinstance(path_of_file, list):
            path_str = str(path_of_file[0])
        else:
            path_str = str(path_of_file)
        
        # More efficient: find last occurrence and extract
        if '.h5' in path_str:
            path_str = path_str.rsplit('.h5', 1)[0]
        if '-' in path_str:
            timestep_str = path_str.rsplit('-', 1)[-1].split('_')[0]
        else:
            timestep_str = path_str.split('_')[0]
        
        return int(timestep_str)
    
    def helper_function(self, all_files):
        # Use set for O(1) lookup instead of list
        timesteps_set = set()
        for file in all_files:
            timestep = self.get_timestep_of_simulation(file)
            timesteps_set.add(timestep)
        timesteps = sorted(timesteps_set)
        # Format once instead of converting to string twice
        str_timesteps = [f"{timestep:05d}" for timestep in timesteps]
        return str_timesteps

    def convert_fields(self):
        """
        method to convert field data from 2D to 3D.
        the method calls the convert_and_write_hdf5_file_fields method from the utilities class in the utils.py file
        it also passes the field, the mode and the domain as arguments to the method from the command line flags


        Parameters:
        -----------
        None

        Returns:
        --------
        None
        """

        # Absolute path to simulation folder
        folder = self.file_path
        uts = self.utils
        uts.set_step_z(self.step_z)
        uts.set_step_r(self.step_r)

        # Define the field to convert ('e1_cyl_m', 'b3_cyl_m,savg', etc.)
        key = f'{self.field}_cyl_m'
        s_avg = ''
        if self.spatial_average :
            key += '-savg'
            s_avg = '-savg'

        if self.lineout:
            key += '-line'
            s_avg = '-line'

        # Define the mode to transform
        mode = self.mode

        # Define the domain size
        if self.domain is None:
            uts.set_domain_size(self.domain)
        elif self.domain is not None:
            domain = {}
            domain['r'] = [self.domain[0],self.domain[1]]
            domain['z'] = [self.domain[2],self.domain[-1]]
            uts.set_domain_size(domain)

        file_extension = 'h5'
        all_files = uts.find_files_by_extension(folder, file_extension)
        filtered_files = uts.filter_files(all_files, self.species, self.data_type)
        print(f"filtered files: {filtered_files}")
        filtered_files = sorted(filtered_files)

        # --- Filter for the specific FIELD this time (e.g., e1, e2, b3) ---
        # Pre-compute field string to avoid repeated f-string formatting
        field_pattern = f'{self.field}_cyl_m'
        files_field = [file for file in filtered_files if field_pattern in file]

        if not files_field:
            print(f"No matching files found for field: {self.field}")
            return
        
        # Now it's safe to access files_field[0]
        source_directory = os.path.dirname(files_field[0]) + '/' #More robust way to create the directory

        # Use os.path.join for more efficient path construction
        source_directory = os.path.dirname(files_field[0])
        if source_directory and not source_directory.endswith('/'):
            source_directory += '/'

        #make the list of fields corresponding to the mode and real or imaginary:
        
        # Determine output directory: if output_path is default ('.'), use source directory
        if self.output_path == '.' or self.output_path == os.getcwd():
            # Create subdirectory in FLD directory: e.g., MS/FLD/e2_3D/
            # Find the FLD directory level
            source_dir = os.path.dirname(files_field[0])
            parts = source_dir.split(os.sep)
            # Find FLD directory and create field_3D subdirectory there
            fld_idx = -1
            for i, part in enumerate(parts):
                if part == 'FLD':
                    fld_idx = i
                    break
            if fld_idx >= 0:
                fld_dir = os.sep.join(parts[:fld_idx+1])
                dest_directory = os.path.join(fld_dir, f'{self.field}_3D')
            else:
                # Fallback: use parent directory
                dest_directory = os.path.join(os.path.dirname(source_dir), f'{self.field}_3D')
        else:
            # Use specified output path
            dest_directory = os.path.join(self.output_path, f'{self.data_type}_{self.field}_to_3D')
        
        if os.path.exists(dest_directory):
            shutil.rmtree(dest_directory)
        os.makedirs(dest_directory, exist_ok=True)

        print(f"Destination directory: {dest_directory}")


        # Loop on all files

        init_time = timeit.default_timer()
        conversion_map = {
            f'e1_cyl_m{s_avg}': '1',
            f'e2_cyl_m{s_avg}': '3',
            f'e3_cyl_m{s_avg}': '4',
            f'b1_cyl_m{s_avg}': '1',
            f'b2_cyl_m{s_avg}': '3',
            f'b3_cyl_m{s_avg}': '4',
            }
        

        timestep = self.timestep
        if timestep is not None and len(timestep) == 1:
            timestep = timestep[0]
        elif timestep is not None and len(timestep) == 2:
            range_timestep = np.arange(timestep[0],timestep[1],1)        
        

        def helper_function(file,dest_directory,conversion_map,mode,i_file,s_avg):

            if (key in conversion_map) and (key.split('_')[0][0] == 'e'):
                uts.convert_and_write_hdf5_file_fields(file,
                                                           dest_directory,
                                                           f'e1_cyl_m{s_avg}', 
                                                           f'e2_cyl_m{s_avg}', 
                                                           f'e3_cyl_m{s_avg}', 
                                                           i_file, 
                                                           mode, 
                                                           conversion_map[key])

            elif (key in conversion_map) and (key.split('_')[0][0] == 'b'):
                uts.convert_and_write_hdf5_file_fields(file,
                                                           dest_directory,
                                                           f'b1_cyl_m{s_avg}', 
                                                           f'b2_cyl_m{s_avg}', 
                                                           f'b3_cyl_m{s_avg}', 
                                                           i_file, 
                                                           mode, 
                                                           conversion_map[key])
            return
        
        if self.field[0] == 'e':
            field_type = 'e'
        else:
            field_type = 'b'
            
        
        def helper_function2(all_files):
            # Use set for O(1) lookup instead of list
            dumps_set = set()
            for file in all_files:
                dump = self.get_timestep_of_simulation(file)
                if dump is not None:
                    dumps_set.add(dump)
            # Format once instead of converting to string twice
            return [f"{dump:06d}" for dump in sorted(dumps_set)]


        all_dumps = helper_function2(filtered_files)
        full_dictionary = {}

        
        for dump in all_dumps:
            #we get the files in this directory
            curr_timestep_files = [file for file in filtered_files if f'{dump}' in file]
            
            
            mode_dict = {}
            for mode in range(self.mode + 1):
                file_dictionary = {}
                comp_dict = {}
                if mode == 0:
                
                    file_dictionary[f'mode_{mode}_{field_type}1_re'] = [file for file in curr_timestep_files
                                                                if f'MODE-{mode}' in file and f'{field_type}1_cyl_m' in file]
                    file_dictionary[f'mode_{mode}_{field_type}2_re'] = [file for file in curr_timestep_files
                                                                if f'MODE-{mode}' in file and f'{field_type}2_cyl_m' in file]
                    file_dictionary[f'mode_{mode}_{field_type}3_re'] = [file for file in curr_timestep_files
                                                                if f'MODE-{mode}' in file and f'{field_type}3_cyl_m' in file]

                    comp_dict['re'] = [file_dictionary[key] for key in file_dictionary.keys() if 're' in key]
                    
                if mode >= 1:
                    # Separate real and imaginary files first
                    re_files = [file for file in curr_timestep_files if f'MODE-{mode}-RE' in file]
                    im_files = [file for file in curr_timestep_files if f'MODE-{mode}-IM' in file]

                    # Assign separately
                    file_dictionary[f'mode_{mode}_{field_type}1_re'] = [file for file in re_files if f'{field_type}1_cyl_m' in file]
                    file_dictionary[f'mode_{mode}_{field_type}2_re'] = [file for file in re_files if f'{field_type}2_cyl_m' in file]
                    file_dictionary[f'mode_{mode}_{field_type}3_re'] = [file for file in re_files if f'{field_type}3_cyl_m' in file]

                    file_dictionary[f'mode_{mode}_{field_type}1_im'] = [file for file in im_files if f'{field_type}1_cyl_m' in file]
                    file_dictionary[f'mode_{mode}_{field_type}2_im'] = [file for file in im_files if f'{field_type}2_cyl_m' in file]
                    file_dictionary[f'mode_{mode}_{field_type}3_im'] = [file for file in im_files if f'{field_type}3_cyl_m' in file]

                    
                    comp_dict['im'] = [file_dictionary[key] for key in file_dictionary.keys() if 'im' in key]
                    comp_dict['re'] = [file_dictionary[key] for key in file_dictionary.keys() if 're' in key]  
                mode_dict[f'mode_{mode}'] = comp_dict
            full_dictionary[f'{dump}'] = mode_dict
         
        # Helper function to check if a timestep should be processed
        def should_process_timestep_fields(dump_val, ts, ts_range):
            if ts is None:
                return True
            elif isinstance(ts, (int, float)):
                return float(dump_val) == ts
            elif isinstance(ts_range, np.ndarray):
                return float(dump_val) in ts_range
            return False

        # Prepare list of dumps to process
        dumps_to_process = []
        for i_file, dump in enumerate(full_dictionary):
            if should_process_timestep_fields(dump, timestep, range_timestep if 'range_timestep' in locals() else None):
                dumps_to_process.append((i_file, dump))

        # Process dumps (parallel or sequential)
        if self.parallel and len(dumps_to_process) > 1:
            print(f'Processing {len(dumps_to_process)} timesteps in parallel with {self.n_workers} workers...')
            args_list = [
                (dump, full_dictionary[dump], dest_directory, conversion_map, mode, i_file, s_avg,
                 self.step_z, self.step_r, self.domain, field_type, key)
                for i_file, dump in dumps_to_process
            ]
            
            with mp.Pool(processes=self.n_workers) as pool:
                results = pool.map(process_single_dump_field, args_list)
            
            # Check for errors
            errors = [r for r in results if not r[1]]
            if errors:
                print(f'Warning: {len(errors)} timestep(s) had errors:')
                for dump, success, error in errors:
                    print(f'  Timestep {dump}: {error}')
        else:
            # Sequential processing (original code)
            for i_file, dump in dumps_to_process:
                print('Converting dump {:d} '.format(i_file))
                helper_function(full_dictionary[dump], dest_directory, conversion_map, mode, i_file, s_avg)

        end_time = timeit.default_timer() - init_time
        print('job finished, total time: ', end_time, 's')

    def convert_raw(self):
            """
            method to convert charge data from 2D to 3D.
            The method calls the convert_and_write_hdf5_file_raws method from the utilities class in the utils.py file
            Parameters:
            -----------
            None

            Returns:
            --------
            None
            """
            
            uts = self.utils
            uts.set_step_z(self.step_z)
            uts.set_step_r(self.step_r)
            folder = self.file_path
            spc = self.species


            # Looking for all files
            file_extension = 'h5'
            all_files = uts.find_files_by_extension(folder,file_extension)
            filtered_files = sorted(uts.filter_files(all_files,self.species,self.data_type))
            #remove to_3D files
            filtered_files = [file for file in filtered_files if 'to_3D' not in file]
            source_directory = os.path.dirname(filtered_files[0])
            if source_directory and not source_directory.endswith('/'):
                source_directory += '/'

            """
            source_directory = folder + 'MS/RAW/' + spc + '/'
            files = sorted( os.listdir(source_directory) )
            N_files = len(files)
            """

            # Determine output directory: if output_path is default ('.'), use source directory
            if self.output_path == '.' or self.output_path == os.getcwd():
                # Create subdirectory in source directory: e.g., MS/RAW/electrons_3D/
                dest_directory = os.path.join(os.path.dirname(source_directory.rstrip('/')), f'{self.species}_3D')
            else:
                # Use specified output path
                dest_directory = os.path.join(self.output_path, f'{self.data_type}_{self.species}_to_3D')
            
            if os.path.exists(dest_directory):
                shutil.rmtree(dest_directory)
            os.makedirs(dest_directory, exist_ok=True)

            print(f"Destination directory: {dest_directory}")
            
            
            if self.domain is None:
                uts.set_domain_size(self.domain)

            elif self.domain is not None:
                domain = {}
                domain['r'] = self.domain[:2]
                domain['z'] = self.domain[2:]
                uts.set_domain_size(domain)
            
            # Loop on all files
            timestep = self.timestep
            if timestep is not None and len(timestep) == 1:
                timestep = timestep[0]
            elif timestep is not None and len(timestep) == 2:
                range_timestep = np.arange(timestep[0],timestep[1],1)
            
            
        
            
            # Helper function to check if a file should be processed
            def should_process_file_raw(file_path, ts, ts_range):
                file_ts = self.get_timestep_of_simulation(file_path)
                if ts is None:
                    return True
                elif isinstance(ts, (int, float)):
                    return file_ts == ts
                elif isinstance(ts_range, np.ndarray):
                    return file_ts in ts_range
                return False

            # Prepare list of files to process
            files_to_process = []
            for i_file, file in enumerate(filtered_files):
                if should_process_file_raw(file, timestep, range_timestep if 'range_timestep' in locals() else None):
                    files_to_process.append((i_file, file))

            # Process files (parallel or sequential)
            if self.parallel and len(files_to_process) > 1:
                print(f'Processing {len(files_to_process)} files in parallel with {self.n_workers} workers...')
                args_list = [
                    (file_path, dest_directory, spc, i_file, self.step_z, self.step_r, self.domain)
                    for i_file, file_path in files_to_process
                ]
                
                with mp.Pool(processes=self.n_workers) as pool:
                    results = pool.map(process_single_dump_raw, args_list)
                
                # Check for errors
                errors = [r for r in results if not r[1]]
                if errors:
                    print(f'Warning: {len(errors)} file(s) had errors:')
                    for i_file, success, error in errors:
                        print(f'  File {i_file}: {error}')
            else:
                # Sequential processing (original code)
                for i_file, file in files_to_process:
                    print('Converting dump {:d} '.format(i_file))
                    uts.convert_and_write_hdf5_file_raws(file, dest_directory, spc, i_file, 
                                                        if_conv_ene_to_gev=True, if_reduce=False, step=1)

            print('job finished')



    def convert_charge(self):
        """
        method to convert raw data from 2D to 3D.
        The method calls the convert_and_write_hdf5_file_densities method from the utilities class in the utils.py file

        Parameters:
        -----------
        None

        Returns:
        --------
        None
        
        """

        # Absolute path to simulation folder
        folder = self.file_path
        uts = self.utils
        uts.set_step_z(self.step_z)
        uts.set_step_r(self.step_r)
        

        # Define the charge to convert
        key = 'charge_cyl_m'
        
        if self.spatial_average :
            key += '-savg'
        
        if self.lineout:
            key += '-line'
        

        # Define the mode to transform
        mode = self.mode
        
        if self.domain is None:
            uts.set_domain_size(self.domain)
        
        if self.domain is not None:
            domain = {}
            domain['r'] = self.domain[:2]
            domain['z'] = self.domain[2:]
            uts.set_domain_size(domain)


        # Looking for all files, This could be a function to avoid code repetition
        # File finding and filtering (using the improved method)
        file_extension = 'h5'
        all_files = uts.find_files_by_extension(folder, file_extension)
        filtered_files = uts.filter_files(all_files, self.species, self.data_type)
        filtered_files = sorted(filtered_files)  # Sorting is good practice


        # You probably don't need this anymore, since filtering is robust
        if filtered_files:  # Avoid IndexError if no files are found
            source_directory = os.path.dirname(filtered_files[0])
            if source_directory and not source_directory.endswith('/'):
                source_directory += '/'
            print(f"Source directory: {source_directory}")

        # ... (Rest of your conversion logic, using filtered_files) ...
        # Make sure you handle the case where filtered_files might be empty!
        if not filtered_files:
            print("No matching files found!")
            return  # Or raise an exception, depending on your needs

        
        """
        source_directory = folder + 'MS/DENSITY/' + self.species + '/MODE-' + str(mode) + '-RE/' + key + '/'
        files = sorted( os.listdir(source_directory))
        files = [f for f in files if not f.startswith('.')]
        N_files = len(files)
        """

        # Determine output directory: if output_path is default ('.'), use source directory
        if self.output_path == '.' or self.output_path == os.getcwd():
            # Create subdirectory in source directory: e.g., MS/DENSITY/electrons_3D/
            # Find the species directory (parent of MODE-X-RE directories)
            source_dir = os.path.dirname(filtered_files[0])
            # Go up to the species level (e.g., from MODE-0-RE/charge_cyl_m/ to electrons/)
            # Keep going up until we're at the species directory level
            parts = source_dir.split(os.sep)
            # Find the index of the species name in the path
            species_idx = -1
            for i, part in enumerate(parts):
                if part == self.species:
                    species_idx = i
                    break
            if species_idx >= 0:
                # Reconstruct path up to species directory
                species_dir = os.sep.join(parts[:species_idx+1])
                dest_directory = os.path.join(species_dir, f'{self.species}_3D')
            else:
                # Fallback: use parent of source directory
                dest_directory = os.path.join(os.path.dirname(source_dir), f'{self.species}_3D')
        else:
            # Use specified output path
            dest_directory = os.path.join(self.output_path, f'{self.data_type}_{self.species}_to_3D')
        
        if os.path.exists(dest_directory):
            shutil.rmtree(dest_directory)
        os.makedirs(dest_directory, exist_ok=True)

        print(f"Destination directory: {dest_directory}")
        


        #check the timestep range. if the timestep is not specified, convert all files, 
        #otherwise convert only the files with the timestep specified or the range of timesteps specified


        # Loop on all files
        timestep = self.timestep
        if timestep is not None and len(timestep) == 1:
            timestep = timestep[0]
        elif timestep is not None and len(timestep) == 2:
            range_timestep = np.arange(timestep[0],timestep[1]+1,1)

        
        def helper_function2(all_files):
            # Use set for O(1) lookup instead of list
            dumps_set = set()
            for file in all_files:
                dump = self.get_timestep_of_simulation(file)
                if dump is not None:
                    dumps_set.add(dump)
            # Format once instead of converting to string twice
            return [f"{dump:06d}" for dump in sorted(dumps_set)]


        
        
        #Here is the good stuff
        all_dumps = helper_function2(filtered_files)
        
        # Auto-detect the maximum mode number from the files
        def detect_max_mode(files):
            max_mode = -1
            for file_path in files:
                # Look for MODE-X patterns in the file path
                mode_matches = re.findall(r'MODE-(\d+)', file_path)
                for mode_str in mode_matches:
                    mode_num = int(mode_str)
                    if mode_num > max_mode:
                        max_mode = mode_num
            return max_mode if max_mode >= 0 else 0  # Default to 0 if no modes found
        
        # Detect max mode from filtered files
        detected_max_mode = detect_max_mode(filtered_files)
        if detected_max_mode > self.mode:
            print(f"Detected mode {detected_max_mode} in files, but user specified mode {self.mode}. Using detected mode {detected_max_mode}.")
            self.mode = detected_max_mode
        elif self.mode == 1 and detected_max_mode == 0:
            print(f"Detected only mode 0 in files, but user specified mode {self.mode}. Using detected mode 0.")
            self.mode = detected_max_mode
        else:
            print(f"Using mode {self.mode} (detected max: {detected_max_mode})")
        
        full_dictionary = {}




        for dump in all_dumps:
            #we get the files in this directory
            curr_timestep_files = [file for file in filtered_files if f'{dump}' in file]
            
            
            mode_dict = {}
            for mode in range(self.mode + 1):
                file_dictionary = {}
                comp_dict = {}
                if mode == 0:

                    file_dictionary[f'mode_{mode}_charge_re'] = [file for file in curr_timestep_files
                                                                if f'MODE-{mode}' in file and f'charge_cyl_m-{self.species}' in file]
                    
                    comp_dict['re'] = [file_dictionary[key] for key in file_dictionary.keys() if 're' in key]
                    
                    
                if mode >= 1:
                    file_dictionary[f'mode_{mode}_charge_re'] = [file for file in curr_timestep_files
                                                                if f'MODE-{mode}' in file and f'charge_cyl_m-{self.species}' in file and 'RE' in file]
                    file_dictionary[f'mode_{mode}_charge_im'] = [file for file in curr_timestep_files
                                                                if f'MODE-{mode}' in file and f'charge_cyl_m-{self.species}' in file and 'IM' in file]

                    
                    comp_dict['im'] = [file_dictionary[key] for key in file_dictionary.keys() if 'im' in key]
                    comp_dict['re'] = [file_dictionary[key] for key in file_dictionary.keys() if 're' in key]  
                mode_dict[f'mode_{mode}'] = comp_dict
            full_dictionary[f'{dump}'] = mode_dict

    

        
        
        #We check if the timestep is None, if it is, we convert all files, otherwise we convert only the files with the timestep specified

        # Helper function to check if a timestep should be processed
        def should_process_timestep(dump_val, ts, ts_range):
            if ts is None:
                return True
            elif isinstance(ts, (int, float)):
                return float(dump_val) == ts
            elif isinstance(ts_range, np.ndarray):
                return float(dump_val) in ts_range
            return False

        # Prepare list of dumps to process
        dumps_to_process = []
        for i_file, dump in enumerate(full_dictionary):
            if should_process_timestep(dump, timestep, range_timestep if 'range_timestep' in locals() else None):
                dumps_to_process.append((i_file, dump))

        # Process dumps (parallel or sequential)
        if self.parallel and len(dumps_to_process) > 1:
            print(f'Processing {len(dumps_to_process)} timesteps in parallel with {self.n_workers} workers...')
            args_list = [
                (dump, full_dictionary[dump], dest_directory, self.species, i_file, mode, 
                 self.step_z, self.step_r, self.domain)
                for i_file, dump in dumps_to_process
            ]
            
            with mp.Pool(processes=self.n_workers) as pool:
                results = pool.map(process_single_dump_charge, args_list)
            
            # Check for errors
            errors = [r for r in results if not r[1]]
            if errors:
                print(f'Warning: {len(errors)} timestep(s) had errors:')
                for dump, success, error in errors:
                    print(f'  Timestep {dump}: {error}')
        else:
            # Sequential processing (original code)
            for i_file, dump in dumps_to_process:
                print('Converting dump {:d} '.format(i_file))
                uts.convert_and_write_hdf5_file_densities(full_dictionary[dump], dest_directory, self.species, i_file, mode)

        print('job finished')
    



def argparser():
    """
    Parser function to specify the input file location and the type of data to be converted.
    The user can also specify the output path, the species, the field, the mode, the domain and the timestep as flags in the command line.

    Flags
        ----------
        file_path : str
            Location of the simulation files

        data_type : str
            Type of data to be converted. The options are charge, field and raw.    
        
        output_path : str
            Location to output the converted files
        
        species : str
            The species to be converted. The default is electrons.
        
        field : str
            The field to be converted. The default is e2.
        
        mode : int
            The mode to be converted. The default is 1.

        domain : list[float]
            The domain to be converted. The default is None.

        timestep : int
            The timestep to be converted. The default is None.

        Returns
        -------
        None
    """
    


    parser = argparse.ArgumentParser(description='Decide which 3D process to do which files')
    parser.add_argument('file_path', metavar='INPUT', type=str, help='input files')
    parser.add_argument('data_type', metavar='DATA_TYPE', type=str, help='type of data to be made 3D')
    parser.add_argument('-o', '--output_path', metavar='OUTPUT', type=str, default='.',
                        help='output file')
    parser.add_argument('-s', '--species', metavar='SPECIES', type=str, default='')
    parser.add_argument('-f', '--field', metavar='FIELD', type=str, default='e2',help='field to be converted')
    parser.add_argument('-m', '--mode', metavar='MODE', type=int, default=1, help='mode to be converted')
    parser.add_argument('-d', '--domain', metavar='DOMAIN', type=float, nargs=4, default=None,help='domain to be converted')
    parser.add_argument('-t', '--timestep', metavar='TIME', type=float,nargs= '+',default=None,help='timestep to be converted')
    parser.add_argument('-savg', '--spatial_average', metavar='SPATIAL_AVERAGE', type=bool, default=False)
    parser.add_argument('-line', '--lineout', metavar='LINEOUT', type=bool, default=False)
    parser.add_argument('-s_z', '--step_z', metavar='STEP_Z',type=int,default=2, help='step in z direction')
    parser.add_argument('-s_r', '--step_r', metavar='STEP_r',type=int,default=2, help='step in in r direction')
    parser.add_argument('-p', '--parallel', action='store_true', help='enable parallel processing of timesteps')
    parser.add_argument('-n', '--n_workers', metavar='N_WORKERS', type=int, default=None, 
                        help='number of parallel workers (default: min(CPU_count, 4))')
    return parser.parse_args()


def main():
    argparse = argparser()
    process_data = ProcessData(argparse.file_path,
                                argparse.data_type,
                                argparse.output_path,
                                argparse.species,
                                argparse.field,
                                argparse.mode,
                                argparse.domain,
                                argparse.timestep,
                                argparse.spatial_average,
                                argparse.lineout,
                                argparse.step_z,
                                argparse.step_r,
                                argparse.parallel,
                                argparse.n_workers)

    
    if argparse.data_type == 'charge':
        print('processing charge data ...')
        process_data.convert_charge()

    elif argparse.data_type == 'field':
        print('processing field data ...')
        process_data.convert_fields()

    elif argparse.data_type == 'raw':
        print('processing raw data ...')
        process_data.convert_raw()
    
    else:
        print('data type not recognized, please choose from the following: charge, field, raw')
        sys.exit(1)

if __name__ == '__main__':
    main()