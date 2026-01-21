# Quick Start: Adding Parallel Processing

## Option 1: Parallel Timesteps (Easiest, Highest Impact)

Add this to `transform_data.py`:

```python
import multiprocessing as mp
from functools import partial

def process_single_dump(args):
    """Process a single dump - for parallel execution"""
    dump, file_dict, dest_directory, species, i_file, mode, utils_instance = args
    print(f'Converting dump {i_file} (timestep {dump})')
    utils_instance.convert_and_write_hdf5_file_densities(
        file_dict, dest_directory, species, i_file, mode
    )
    return dump

# In convert_charge(), replace the sequential loop with:
if self.parallel and len(full_dictionary) > 1:
    n_workers = min(mp.cpu_count(), self.n_workers or 4)
    args_list = [
        (dump, full_dictionary[dump], dest_directory, self.species, i_file, mode, uts)
        for i_file, dump in enumerate(full_dictionary)
        if self.should_process_timestep(dump, timestep, range_timestep)
    ]
    
    with mp.Pool(processes=n_workers) as pool:
        pool.map(process_single_dump, args_list)
else:
    # Original sequential code
    for i_file, dump in enumerate(full_dictionary):
        # ... existing code ...
```

## Option 2: Parallel Z-Slices (Moderate Complexity)

Modify `utils.pyx` to process slices in parallel using ThreadPoolExecutor:

```python
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp

def process_slices_parallel(self, data, r, R, theta, mode_number, nz):
    """Process z-slices in parallel"""
    def process_slice(i):
        interpolated = self.compute_density_for_slice(i, data, r, R, theta, mode_number)
        return i, np.nan_to_num(interpolated, copy=False)
    
    # Process in parallel
    with ThreadPoolExecutor(max_workers=min(mp.cpu_count(), nz)) as executor:
        results = list(executor.map(process_slice, range(nz)))
    
    # Write sequentially to avoid HDF5 conflicts
    with h5py.File(self.target_path, 'a') as f:
        dataset = f['/dataset']
        for i, output_slice in sorted(results):
            dataset[:, :, i] = output_slice
```

## Option 3: Optimize Interpolation (Easy, Medium Impact)

Replace repeated `np.interp` calls with pre-computed interpolators:

```python
# Pre-compute interpolators once
interpolators = {}
for mode in range(mode_number + 1):
    for part in ['re', 'im']:
        key = f'mode_{mode}_{part}_charge'
        if key in data:
            interpolators[key] = interp1d(r, data[key], 
                                         kind='linear', 
                                         bounds_error=False, 
                                         fill_value=0.0)

# Then use in compute_density_for_slice:
charge_re_interp = interpolators[key_re](R)
```

## Expected Performance Gains

- **Parallel timesteps**: 2-4x speedup (4 cores) to 4-8x (8+ cores)
- **Parallel slices**: 2-3x speedup depending on slice count
- **Optimized interpolation**: 1.5-2x speedup
- **Combined**: 5-15x total speedup possible

## Memory Considerations

- Each parallel worker uses its own memory
- For large datasets, limit workers: `n_workers = min(4, mp.cpu_count())`
- Consider processing timesteps sequentially but slices in parallel

## Testing

Test with a small dataset first:
```bash
python transform_data.py examples/beam.dev/MS/ charge -s electrons -t 0 1
```

Then scale up to full dataset.
