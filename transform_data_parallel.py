"""
Parallel processing version of transform_data.py
This demonstrates how to add multiprocessing for timesteps and slices
"""

import multiprocessing as mp
from functools import partial
import numpy as np
from transform_data import ProcessData
import utils

def process_single_timestep(args):
    """
    Process a single timestep - designed to be called in parallel
    """
    dump, file_dict, dest_directory, species, i_file, mode, utils_instance = args
    
    print(f'Converting dump {i_file} (timestep {dump})')
    try:
        utils_instance.convert_and_write_hdf5_file_densities(
            file_dict, dest_directory, species, i_file, mode
        )
        return (dump, True, None)
    except Exception as e:
        return (dump, False, str(e))


def process_slices_parallel(utils_instance, data, r, R, theta, mode_number, nz, output_path):
    """
    Process z-slices in parallel
    Note: This requires careful handling of HDF5 file access
    """
    from concurrent.futures import ThreadPoolExecutor
    import h5py
    
    def process_single_slice(i):
        """Process a single z-slice"""
        interpolated_densities = utils_instance.compute_density_for_slice(
            i, data, r, R, theta, mode_number
        )
        output_slice = np.nan_to_num(interpolated_densities, copy=False)
        return i, output_slice
    
    # Process slices in parallel
    with ThreadPoolExecutor(max_workers=mp.cpu_count()) as executor:
        results = list(executor.map(process_single_slice, range(nz)))
    
    # Write results sequentially to avoid HDF5 conflicts
    with h5py.File(output_path, 'a') as f:
        dataset = f['/dataset']
        for i, output_slice in results:
            dataset[:, :, i] = output_slice


class ProcessDataParallel(ProcessData):
    """
    Parallel version of ProcessData that processes timesteps in parallel
    """
    
    def convert_charge_parallel(self, n_workers=None):
        """
        Parallel version of convert_charge that processes multiple timesteps simultaneously
        """
        if n_workers is None:
            n_workers = min(mp.cpu_count(), 4)  # Limit to 4 to avoid memory issues
        
        # Get all the setup from parent class
        folder = self.file_path
        uts = self.utils
        uts.set_step_z(self.step_z)
        uts.set_step_r(self.step_r)
        
        # ... (same setup code as convert_charge) ...
        
        # Prepare arguments for parallel processing
        args_list = []
        for i_file, dump in enumerate(full_dictionary):
            if timestep is None or self.should_process_timestep(dump, timestep, range_timestep):
                args_list.append((
                    dump,
                    full_dictionary[dump],
                    dest_directory,
                    self.species,
                    i_file,
                    mode,
                    uts
                ))
        
        # Process in parallel
        print(f'Processing {len(args_list)} timesteps with {n_workers} workers...')
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(process_single_timestep, args_list)
        
        # Check for errors
        for dump, success, error in results:
            if not success:
                print(f'Error processing timestep {dump}: {error}')
        
        print('job finished')


# Usage example:
# process_data = ProcessDataParallel(...)
# process_data.convert_charge_parallel(n_workers=4)
