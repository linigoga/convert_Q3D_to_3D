# Parallel Processing Usage Guide

## Overview

Parallel timestep processing has been implemented to speed up conversion of multiple timesteps. Each timestep/dump is processed independently in parallel, providing significant speedup for datasets with many timesteps.

## Usage

### Basic Usage (Sequential - Default)

```bash
python transform_data.py examples/beam.dev/MS/ charge -s electrons
```

### Enable Parallel Processing

Add the `-p` or `--parallel` flag:

```bash
python transform_data.py examples/beam.dev/MS/ charge -s electrons -p
```

### Specify Number of Workers

Control the number of parallel workers with `-n` or `--n_workers`:

```bash
python transform_data.py examples/beam.dev/MS/ charge -s electrons -p -n 8
```

**Default**: `min(CPU_count, 4)` - limits to 4 workers by default to avoid memory issues

### Examples

1. **Charge conversion with parallel processing:**
   ```bash
   python transform_data.py examples/beam.dev/MS/ charge -s electrons -p
   ```

2. **Field conversion with 8 workers:**
   ```bash
   python transform_data.py examples/beam.dev/MS/ field -f e2 -m 1 -p -n 8
   ```

3. **Raw data conversion with parallel processing:**
   ```bash
   python transform_data.py examples/beam.dev/MS/ raw -s electrons -p
   ```

4. **Specific timestep range with parallel processing:**
   ```bash
   python transform_data.py examples/beam.dev/MS/ charge -s electrons -t 0 10 -p
   ```

## Performance Expectations

- **2-4 CPU cores**: 2-3x speedup
- **4-8 CPU cores**: 4-6x speedup  
- **8+ CPU cores**: 6-8x speedup (limited by I/O)

**Note**: Speedup depends on:
- Number of timesteps to process
- Dataset size (larger datasets benefit more)
- Disk I/O speed
- Available memory

## Memory Considerations

- Each worker process uses its own memory space
- For large datasets, limit workers: `-n 2` or `-n 4`
- Monitor memory usage: `htop` or `top`
- If you run out of memory, reduce `-n` value

## When to Use Parallel Processing

✅ **Use parallel processing when:**
- Processing multiple timesteps (> 2)
- Have sufficient RAM (8GB+ recommended)
- CPU has multiple cores
- Dataset files are not extremely large

❌ **Avoid parallel processing when:**
- Processing single timestep
- Limited RAM (< 4GB)
- Very large individual files (> 10GB each)
- On systems with limited CPU cores

## Troubleshooting

### Out of Memory Errors

Reduce the number of workers:
```bash
python transform_data.py ... -p -n 2
```

### HDF5 File Locking Issues

If you see file locking errors, the code handles this automatically by writing sequentially even when processing in parallel. If issues persist, reduce workers or process sequentially.

### Progress Monitoring

The code prints progress for each worker:
```
Processing 10 timesteps in parallel with 4 workers...
Converting dump 0 (timestep 000000)
Converting dump 1 (timestep 000001)
...
```

## Implementation Details

- Uses Python's `multiprocessing.Pool` for parallel execution
- Each timestep is processed in a separate process
- Processes are independent (no shared state)
- Error handling: errors in one timestep don't stop others
- Automatic worker cleanup after completion

## Backward Compatibility

The code is fully backward compatible:
- Without `-p` flag: runs sequentially (original behavior)
- With `-p` flag: runs in parallel
- All existing command-line arguments work the same way
