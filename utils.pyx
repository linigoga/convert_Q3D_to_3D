import h5py
import shutil
import os
import numpy as np  # For general numpy functions
import glob # For finding files in a directory

from scipy.interpolate import interp1d  # Regular import
from scipy.ndimage import map_coordinates  # Regular import
from scipy.interpolate import CubicSpline  # Regular import

cimport numpy as np  # For Cython specific features like typed memoryvie


cdef class utilities:
    cdef dict domain_size
    cdef str target_path
    cdef int step_z
    cdef int step_r
    cdef tuple files_to_read


    def __init__(self) -> None:
        self.domain_size = {}
        self.target_path = None
        

    cpdef void set_domain_size(self,dict domain_size):
        self.domain_size = domain_size
        return
    
    cpdef void set_target_path(self,str target_path):
        self.target_path = target_path
        return

    cpdef void set_step_z(self,int step_z):
        self.step_z = step_z
        return

    cpdef void set_step_r(self,int step_r):
        self.step_r = step_r
        return
    
    cpdef list find_files_by_extension(self, str target_path, str file_extension):
            
        file_paths = glob.glob(os.path.join(target_path, '**', f'*.{file_extension}'), recursive=True)
        return file_paths


    cpdef list filter_files(self, list file_paths, str species, str data_type):
        """
        Filters files based on species and data type, allowing for flexible directory structures.
        """

        if data_type == 'raw':
            data_type_keyword = 'RAW'
        elif data_type == 'field':
            data_type_keyword = 'FLD'
        elif data_type == 'density' or data_type == 'charge':
            data_type_keyword = 'DENSITY'
        else:
            raise Exception('data_type not recognized')  # Use raise for exceptions

        # First filter: Check if data_type_keyword is in the path

        print("data_type_keyword",data_type_keyword)
        filtered_files = [file for file in file_paths if data_type_keyword in file]
        
        if species != '':
            print('species' + species)
            filtered_files = [file for file in filtered_files if species in file]

        
        # Third filter: remove the 3D converted files
        filtered_files = [file for file in filtered_files if 'to_3D' not in file]

        return filtered_files


    cpdef tuple read_file_2d(self, str path):

        """
        Reads a 2d osiris file and returns the data and the bounds of the grid

        Parameters:
        ---------
        input:

        path: string
        outputs:
        - min_z, max_z, nz: min bound, max bound, number of cells on axis z
        - min_r, max_r, nr: min bound, max bound, number of cells on axis r
        - fdata: the 2D array containing the data, size (nr,nz)
        """
        cdef double min_z, max_z, min_r, max_r
        cdef int nz, nr

        # Open the file
        with h5py.File(path, 'r') as f:

            # Get keys of the file
            f_keys = list(f.keys())
            data = f[f_keys[-1]]
            nr, nz = data.shape[0], data.shape[1]

            # Read the array and store it
            fdata = np.zeros((nr, nz), dtype=np.float64)
            fdata_view = fdata  # <- Assign after creating fdata
            data.read_direct(fdata, np.s_[0:nr, 0:nz], np.s_[0:nr, 0:nz])


            # Get the mins and maxs bounds along each direction
            min_z = f['AXIS/AXIS1'][:][0]
            max_z = f['AXIS/AXIS1'][:][1]
            min_r = f['AXIS/AXIS2'][:][0]
            max_r = f['AXIS/AXIS2'][:][1]

        return min_z, max_z, nz, min_r, max_r, nr, fdata



    cpdef void create_osiris_h5(self, str path):
        #cdef h5py.File f
        #cdef h5py.Group grp
        #cdef h5py.Dataset dset
        
        """
        Creates a template 3d osiris file with default parameters
        input :
        - path : absolute path to the file
        outputs :
        - none
        """  
        
        f = h5py.File(path, "w")
        f.attrs['ITER'] = [np.int32(0)]
        f.attrs['LABEL'] = np.array([b'\\rho '], dtype=h5py.string_dtype('ascii', 256))
        f.attrs['NAME'] = np.array([b'charge '], dtype=h5py.string_dtype('ascii', 256))
        f.attrs['TIME'] = [float(0)]
        f.attrs['TIME UNITS'] = np.array([b'1 / \\omega_p '], dtype=h5py.string_dtype('ascii', 256))
        f.attrs['TYPE'] = np.array([b'grid'], dtype=h5py.string_dtype('ascii', 4))
        f.attrs['UNITS'] = np.array([b'e \\omega_p^3/ c^3 '], dtype=h5py.string_dtype('ascii', 256))
        grp = f.create_group("AXIS")
        dset = grp.create_dataset("AXIS1",  data=[0., 1.], dtype='f8')
        dset.attrs['LONG_NAME'] = np.array([b'x_1 '], dtype=h5py.string_dtype('ascii', 256))
        dset.attrs['NAME'] = np.array([b'x1 '], dtype=h5py.string_dtype('ascii', 256))
        dset.attrs['TYPE'] = np.array([b'linear'], dtype=h5py.string_dtype('ascii', 6))
        dset.attrs['UNITS'] = np.array([b'c / \\omega_p '], dtype=h5py.string_dtype('ascii', 256))
        dset = grp.create_dataset("AXIS2",  data=[0., 1.], dtype='f8')
        dset.attrs['LONG_NAME'] = np.array([b'x_2 '], dtype=h5py.string_dtype('ascii', 256))
        dset.attrs['NAME'] = np.array([b'x2 '], dtype=h5py.string_dtype('ascii', 256))
        dset.attrs['TYPE'] = np.array([b'linear'], dtype=h5py.string_dtype('ascii', 6))
        dset.attrs['UNITS'] = np.array([b'c / \\omega_p '], dtype=h5py.string_dtype('ascii', 256))
        dset = grp.create_dataset("AXIS3",  data=[0., 1.], dtype='f8')
        dset.attrs['LONG_NAME'] = np.array([b'x_3 '], dtype=h5py.string_dtype('ascii', 256))
        dset.attrs['NAME'] = np.array([b'x3 '], dtype=h5py.string_dtype('ascii', 256))
        dset.attrs['TYPE'] = np.array([b'linear'], dtype=h5py.string_dtype('ascii', 6))
        dset.attrs['UNITS'] = np.array([b'c / \\omega_p '], dtype=h5py.string_dtype('ascii', 256))
        grp = f.create_group("SIMULATION")
        grp.attrs['COMPILE_TIME'] = np.array([b'Oct 7 2021 16:01:45 '], dtype=h5py.string_dtype('ascii', 1024))
        grp.attrs['DT'] = [float(0.02)]
        grp.attrs['GIT_VERSION'] = np.array([b'4.4.4-212-g8a4effe '], dtype=h5py.string_dtype('ascii', 1024))
        grp.attrs['INPUT_FILE'] = np.array([b'os-stdin '], dtype=h5py.string_dtype('ascii', 1024))
        grp.attrs['INPUT_FILE_CRC32'] = [float(147908860)]
        grp.attrs['MOVE C'] = [np.int32(0),np.int32(0),np.int32(0)]
        grp.attrs['NDIMS'] = [np.int32(3)]
        grp.attrs['NX'] = [np.int32(100),np.int32(100),np.int32(100)]
        grp.attrs['PAR_NODE_CONF'] = [np.int32(1),np.int32(1),np.int32(1)]
        grp.attrs['PAR_NX_X1'] = [np.int32(100)]
        grp.attrs['PAR_NX_X2'] = [np.int32(100)]
        grp.attrs['PAR_NX_X3'] = [np.int32(100)]
        grp.attrs['PERIODIC'] = [np.int32(0),np.int32(0),np.int32(0)]
        grp.attrs['TIMESTAMP'] = np.array([b''], dtype=h5py.string_dtype('ascii', 1024))
        grp.attrs['XMAX'] = [float(1),float(1),float(1)]
        grp.attrs['XMIN'] = [float(0),float(0),float(0)]
        f.close()

    cpdef void prepare_hdf5_file(self,
                          str path,
                          str target_path,
                          float min_z,
                          float max_z,
                          float max_r, 
                          int nx, 
                          int ny, 
                          int nz):

        '''
        Write the 3D hdf5 file where the cartesian data will be stored
        Parameters :
        ----------
        inputs : 
        path : string
            absolute path to the cylindrical file to convert

        target_path : string
            absolute path to the folder where the converted file will be stored
        min_z : float
            min bound of the cylindrical grid in x1

        max_z : float
            max bound of the cylindrical grid in x1

        max_r : float
            max bound of the cylindrical grid in x4
        
        nx, ny, nz : integers
            number of points for the cartesian grid
        
        outputs :
        ----------
            None
        '''


        cdef str target_file
        
        
        target_file = os.path.join(target_path, path.split('/')[-1].split('.')[0])
        
        dst = target_file + '.h5'
        self.set_target_path(dst)
        
        if os.path.isfile(dst) :
                os.remove(dst)
        
        
        # create template osiris-3d file
        self.create_osiris_h5(dst)

        # We copy the attributes from the file we convert
        src = h5py.File(path, 'r')
        # Read the cartesian 3D template osiris output file
        dst = h5py.File(dst, 'a')

        src_attrs = src.attrs
        dst_attrs = dst.attrs


        
        for name in dst_attrs.keys() :
            dst_attrs.modify(name, src_attrs.__getitem__(name))

        # Change the bounds of each axis
        dst['AXIS/AXIS1'][0:2] = [min_z, max_z]
        dst['AXIS/AXIS2'][0:2] = [-max_r, max_r]
        dst['AXIS/AXIS3'][0:2] = [-max_r, max_r]

        # Create a new dataset with the dimensions
        dst.create_dataset('/dataset', (ny, nx, nz,), dtype=np.double)

        # close files
        src.close()
        dst.close()


        return



    cpdef np.ndarray[np.double_t, ndim=2] polar2cartesian(self,
                    np.ndarray[np.double_t] r, 
                    np.ndarray[np.double_t] t, 
                    np.ndarray[np.double_t, ndim=2] grid, 
                    np.ndarray[np.double_t] x, 
                    np.ndarray[np.double_t] y,
                    int order=3):
        """
        Convert a 2D grid from cylindrical to cartesian coordinates.

        Parameters:
        ----------
        r : 1D array
            Radial coordinate of the cylindrical grid.
        t : 1D array
            Azimuthal coordinate of the cylindrical grid.
        grid : 2D array
            The 2D array containing the data with dimensions (nr, nz).
        x, y : 1D arrays
            x, y coordinates of the cartesian grid.
        order : integer, optional (default: 3)
            Order of the interpolation.

        Returns:
        -------
        cart_data : 2D array
            The cartesian 2D array containing the data with dimensions (nx, ny).
        """
        cdef np.ndarray[np.double_t, ndim=2] cart_data, X, Y, new_r, new_t
        cdef int len_r = len(r)
        cdef int len_t = len(t)

        X, Y = np.meshgrid(x, y)
        new_r = np.sqrt(X**2 + Y**2)
        new_t = np.arctan2(Y, X)
        

        # If r and t are evenly spaced:
        r_spacing = (r[-1] - r[0]) / (len_r - 1)
        t_spacing = (t[-1] - t[0]) / (len_t - 1)

        # Calculate and clip the new indices:
        new_ir = np.clip(((new_r - r[0]) / r_spacing).astype(np.int64), 0, len_r-1)
        new_it = np.clip(((new_t - t[0]) / t_spacing).astype(np.int64), 0, len_t-1)


        # Interpolate the data:
        coordinates = np.vstack([new_ir.ravel(), new_it.ravel()]).reshape(2, *new_ir.shape)

        cart_data = map_coordinates(grid, coordinates, order=order)
        cart_data = cart_data.reshape(np.shape(new_r))

        # Mask data values outside the maximum radial extent of the original cylindrical grid:
        mask = (new_r > r[-1])
        cart_data[mask] = 0.0

        return cart_data


    cpdef void convert_and_write_hdf5_file_raws(self,
                                        str file,
                                        str target_directory, 
                                        str spc, 
                                        int n, 
                                        bint if_conv_ene_to_gev, 
                                        bint if_reduce, 
                                        int step):

        '''
        Convert a RAW file from cylindrical to 3D cartesian
        the file with 3D data is the same as the original + '_converted.h5'
        Parameters :
        ----------
        inputs :
        
        file : string
            absolute path to the cylindrical file to convert
        
        target_path : string
            absolute path to the folder where the converted file will be stored

        spc : string
            species name (electron, ion, etc.)
        
        n : integer
            number of the diagnostic dump we want to convert

        if_conv_ene_to_gev : boolean
            set to True to get the 'ene' field in GeV instead of mc2
        
        if_reduce : boolean
            set to True to select one particle every "step"
        
        step : integer
            to sample particles from the RAW and reduce its size
        
        ----------
        outputs:
            None
        '''
        cdef str src, dst
        #cdef h5py.File f
        cdef char* new_quants[9]
        cdef char* new_labels[9]
        cdef char* new_units[9]


        target_file = os.path.join(target_directory, file.split('/')[-1].split('.')[0])
        
        dst = target_file + '.h5'
        self.set_target_path(dst)

        if os.path.isfile(dst) :
                os.remove(dst)


        # File from quasi-3D simulation
        src = file

        # File where we will convert
        dst = target_directory + file.split('/')[-1]

        target_file = os.path.join(target_directory, file.split('/')[-1].split('.')[0])
        

        dst = target_file + '.h5'
        self.set_target_path(dst)

        print(src,dst)
        # Create the file
        if os.path.isfile(dst):
            os.remove(dst)
        shutil.copyfile(src,dst)
        

        # Open the 3D raw file
        f = h5py.File(dst ,'a')

        # Replace the values of x2 by the values of x4
        try :
            f['x2'][::] = f['x4'][::]
        except ValueError :
            f['x2'][...] = 0.

        # Delete x4
        del f['/x4']

        # Reduce the size of the raw file
        if if_reduce :
            for name in f.__iter__():
                if name != 'SIMULATION':
                    try :
                        tmp = f[name][::step]
                        del f[name]
                        f[name] = tmp
                    except ValueError :
                        f[name][...] = 0.

        # Replace the values of ene by the values of ene in GeV
        if if_conv_ene_to_gev :
            try :
                f['ene'][::] = f['ene'][::] * 0.511e-3
            except ValueError :
                f['ene'][...] = 0.

        # Change the QUANTS
        f.attrs.__delitem__('QUANTS')
        new_quants = [b'x1', b'x2', b'x3', b'p1', b'p2', b'p3', b'q', b'ene', b'tag']
        f.attrs.__setitem__('QUANTS', new_quants)

        # Change the LABELS
        f.attrs.__delitem__('LABELS')
        new_labels = [b'x_1', b'x_2', b'x_3', b'p_1', b'p_2', b'p_3', b'q', b'Ene', b'Tag']
        f.attrs.__setitem__('LABELS', new_labels)

        # Change the UNITS
        f.attrs.__delitem__('UNITS')
        new_units = [b'c/\omega_p', b'c/\omega_p', b'c/\omega_p', b'm_e c', b'm_e c', b'm_e c', b'e', b'm_e c^2', b'']
        f.attrs.__setitem__('UNITS', new_units)

        # Close
        f.close()

        return



    cdef compute_density_for_slice(self, i, data, r, R, t, mode_number):

        '''
        This function computes the charge density for a given slice i.
        It is used in the function convert_and_write_hdf5_file_fields.

        The main tasks include:
        - Interpolating charge data
        - Summing up modes together
        '''

        # Try fetching mode_0 data for the given slice.
        # This acts as our base data upon which other mode contributions are added.
        try:
            charge = data[f'mode_0_re_charge'][:, i]
        except KeyError:
            # If the slice does not exist in our dataset, log an error and return an empty result.
            print('key error, the file is probably empty')
            return np.zeros_like(R)

        # Use np.interp to generate initial interpolated charge
        interpolated_charge = np.interp(R, r, charge)
        
        # Loop through all other modes to sum up their contributions
        for mode in range(1, mode_number + 1): 

            # Fetch and interpolate the real part of the charge for the current mode
            charge_re_interp = np.interp(R, r, data[f'mode_{mode}_re_charge'][:,i])

            # Fetch and interpolate the imaginary part of the charge for the current mode
            charge_im_interp = np.interp(R, r, data[f'mode_{mode}_im_charge'][:,i])

            # Compute cosine and sine values needed to sum real and imaginary parts
            cos_vals = np.cos(mode * t)
            sin_vals = np.sin(mode * t)
            
            # Update the interpolated charge by adding the contribution of the current mode
            interpolated_charge += charge_re_interp * cos_vals + charge_im_interp * sin_vals

        # Return the final sum of interpolated charge
        return interpolated_charge





    cpdef tuple convert_and_write_hdf5_file_densities(self,
                                              object file,
                                              str target_directory, 
                                              str spc,  
                                              int n, 
                                              int mode_number):

        '''
        This function converts osiris output data from cylindrical to 3D cartesian
        the file with 3D data is the same as the original + '_converted.h5'
        
        Parameters :
        ----------
        inputs :
        path : string
            it is the path to a density file output from osiris quasi-3D

        spc : string
            particle species
        
        key : string
            choose "charge_cyl_m"

        n : integer
            diagnostic dump number we want to convert
        
        mode : integer
            0 or 1, to choose which mode we want to convert
        
        ----------
        outputs
            None
        '''

        cdef str path_shortcut
        cdef np.ndarray[np.double_t] axis_z, axis_r, r, t, x, y
        cdef int min_idz, max_idz, min_idr, max_idr, step_z, step_r, nt, nx, ny, nz
        

        #--------------------------------------------
        # Read all files needed
        #--------------------------------------------

        ###Check if the files added here are corresponding to the real and imaginary parts
        data = {}
        for mode, mode_data in file.items():  # mode will be '0', '1', etc.
            for part, paths in mode_data.items():  # part will be 're' or 'im'
                for idx, path in enumerate(paths):  # idx will index over 'path_to_file1', 'path_to_file2', etc.
                    min_z, max_z, nz, min_r, max_r, nr, charge_data = self.read_file_2d(path[0])
                    key = f'{mode}_{part}_charge'
                    data[key] = charge_data  
        
        step_z,step_r = self.step_z, self.step_r
        axis_z = np.linspace(min_z, max_z, nz)
        axis_r = np.linspace(min_r, max_r, nr)

        if self.domain_size is None:
            pass
        else:
            min_z += self.domain_size['z'][0]
            max_z -= self.domain_size['z'][1]
            min_r = self.domain_size['r'][0]
            max_r = self.domain_size['r'][1]

        # indexes in the diagnostic grid
        min_idz = abs(axis_z-min_z).argmin()
        max_idz = abs(axis_z-max_z).argmin()
        min_idr = abs(axis_r-min_r).argmin()
        max_idr = abs(axis_r-max_r).argmin()

        #Have a loop to correct the domain and save on site
        for mode, mode_data in file.items():  # mode will be '0', '1', etc.
            for part, paths in mode_data.items():  # part will be 're' or 'im'
                for idx, path in enumerate(paths):  # idx will index over 'path_to_file1',
                    key = f"{mode}_{part}_charge"
                    data[key] = data[key][min_idr:max_idr:step_r, min_idz:max_idz:step_z]

        

        sample_key = list(data.keys())[0]

        nr, nz = data[sample_key].shape
        
        #--------------------------------------------
        # Define the cylindrical grid
        #--------------------------------------------

        min_t, max_t, nt = -np.pi, np.pi, 8

        r = np.linspace(min_r, max_r, nr)
        t = np.linspace(min_t, max_t, nt)

        #--------------------------------------------
        # Define the cartesian grid
        #--------------------------------------------

        nx = nr
        ny = nr

        x = np.linspace(-max_r, max_r, nx)
        y = np.linspace(-max_r, max_r, ny)


        #print(file[f'mode_{mode_number}']['re'])
        path_shortcut = file[f'mode_{mode_number}']['re'][0][0]

        
        self.prepare_hdf5_file(path_shortcut,target_directory ,min_z, max_z, max_r, nx, ny, nz) # <- check this and how it could be ammended
        
        
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        #print(np.arctan2(Y,X))
        theta = np.arctan2(Y,X)


        #--------------------------------------------
        # Fill the 3D cartesian hdf5 file slice by slice
        #--------------------------------------------

        # Open this file
        f = h5py.File(self.target_path, 'a')

        # Loop on all the z slices
        for i in range(nz):
            
            # Use cylindrical symmetry to have the transverse slice
            # not that only mode 0 is used here
            interpolated_densities = self.compute_density_for_slice(i, data, r, R, theta, mode_number)
            
            
            f['/dataset'][:,:,i] = np.nan_to_num(interpolated_densities)


        f.close()

        return 

    cdef compute_fields_for_slice(self, i, data, r, R, t, mode_number):

        '''
        This function computes the fields for a given slice i.
        It is used in the function convert_and_write_hdf5_file_fields.
        
        The main tasks include:
        - Interpolating fields
        - Summing up modes together
        '''

        # Try fetching mode_0 data for the given slice.
        # This acts as our base data upon which other mode contributions are added.
        try:
            Ez = data[f'mode_0_re_field_1'][:, i]
            Er = data[f'mode_0_re_field_2'][:, i]
            Et = data[f'mode_0_re_field_3'][:, i]
        except IndexError:
            # If the slice does not exist in our dataset, log an error and return an empty result.
            print('index error, the file is probably empty')
            return {}

        # Use np.interp to generate initial interpolated fields
        interpolated_fields = {
            'field_1': np.interp(R, r, Ez),
            'field_2': np.interp(R, r, Er),
            'field_3': np.interp(R, r, Et)
        }

        # Loop through all other modes to sum up their contributions
        for mode in range(1, mode_number + 1): 

            # Fetch and interpolate the real parts of the fields for the current mode
            field_1_re_interp = np.interp(R, r, data[f'mode_{mode}_re_field_1'][:,i])
            field_2_re_interp = np.interp(R, r, data[f'mode_{mode}_re_field_2'][:,i])
            field_3_re_interp = np.interp(R, r, data[f'mode_{mode}_re_field_3'][:,i])

            # Fetch and interpolate the imaginary parts of the fields for the current mode
            field_1_im_interp = np.interp(R, r, data[f'mode_{mode}_im_field_1'][:,i])
            field_2_im_interp = np.interp(R, r, data[f'mode_{mode}_im_field_2'][:,i])
            field_3_im_interp = np.interp(R, r, data[f'mode_{mode}_im_field_3'][:,i])

            # Compute cosine and sine values needed to sum real and imaginary parts
            cos_vals = np.cos(mode * t)
            sin_vals = np.sin(mode * t)
            
            # Update the interpolated fields by adding the contribution of the current mode
            interpolated_fields['field_1'] += field_1_re_interp * cos_vals + field_1_im_interp * sin_vals
            interpolated_fields['field_2'] += field_2_re_interp * cos_vals + field_2_im_interp * sin_vals
            interpolated_fields['field_3'] += field_3_re_interp * cos_vals + field_3_im_interp * sin_vals

        # Return the final sum of interpolated fields
        return interpolated_fields






    cpdef void convert_and_write_hdf5_file_fields(self,
                                           dict file,
                                           str target_directory,
                                           str key_z,
                                           str key_r, 
                                           str key_t, 
                                           int n, 
                                           int mode_number, 
                                           str dir):

        '''
        This function converts osiris output data from cylindrical to 3D cartesian
        
        It works only for fields diagnostics E2, E3, B2, B3 for now
        The file with 3D data is the same as the original + '_converted.h5'

        Parameters :
        ----------
        inputs :
        
        path : list
            list containing the path to the cylindrical file to convert, contains
        
        key_z : string
            choose "e1_cyl_m" or "b1_cyl_m" to convert E1 or B1
        
        key_r : string
            choose "e2_cyl_m", "e3_cyl_m" to convert E2 and E3

        key_t : string
            choose "b2_cyl_m", "b3_cyl_m" to convert B2 and B3

        n : integer
            diagnostic dump number we want to convert   
        
        mode : integer
            0 or 1, to choose which mode we want to convert

        dir : string
            '3' (resp. '4') to specify whether we want E3 (resp. E2)

        ----------
        outputs
            None
        '''

        cdef str path_shortcut
        cdef np.ndarray[np.double_t] axis_z, axis_r, r, t , x, y
        cdef int min_idz, max_idz, min_idr, max_idr, step_z, step_r, nt, nx, ny, nz
        cdef np.ndarray[np.double_t, ndim=2] theta
        


        #--------------------------------------------
        # Read all files needed
        #--------------------------------------------

        #maybe we can simply remove the key, since we are feeding the file? or just keep it, and discriminate later on?

        # the names here could be improved, for instance, intead of Ez, it could be mode_0_rad_re, etc.
        # we can input a the list of fields of this specific timestep?


        data = {}  # This will store the loaded data
        

        for mode, mode_data in file.items():  # mode will be '0', '1', etc.
            for part, paths in mode_data.items():  # part will be 're' or 'im'
                for idx, path in enumerate(paths):  # idx will index over 'path_to_file1', 'path_to_file2', etc.
                    min_z, max_z, nz, min_r, max_r, nr, field_data = self.read_file_2d(path[0])
                    # Generate a key based on mode, part (real/imaginary), and field index (Ez/Er/Et)
                    
                    key = f"{mode}_{part}_field_{idx + 1}"
                    data[key] = field_data  
        
        step_z,step_r = self.step_z, self.step_r
        axis_z = np.linspace(min_z, max_z, nz)
        axis_r = np.linspace(min_r, max_r, nr)
        
        if self.domain_size is None:
            pass
        else:
            min_z += self.domain_size['z'][0]
            max_z -= self.domain_size['z'][1]
            min_r = self.domain_size['r'][0]
            max_r = self.domain_size['r'][1]

        # indexes in the diagnostic grid
        min_idz = abs(axis_z-min_z).argmin()
        max_idz = abs(axis_z-max_z).argmin()
        min_idr = abs(axis_r-min_r).argmin()
        max_idr = abs(axis_r-max_r).argmin()

        #Have a loop to correct the domain and save on site
        for mode, mode_data in file.items():  # mode will be '0', '1', etc.
            for part, paths in mode_data.items():  # part will be 're' or 'im'
                
                for idx, path in enumerate(paths):  # idx will index over 'path_to_file1',
                    key = f"{mode}_{part}_field_{idx + 1}"

                    data[key] = data[key][min_idr:max_idr:step_r, min_idz:max_idz:step_z]
                    print(key)


        
        sample_key = f'mode_0_re_field_1'

        nr, nz = data[sample_key].shape
        
        
        #--------------------------------------------
        # Define the cylindrical grid
        #--------------------------------------------

        min_t, max_t, nt = -np.pi, np.pi, 8    #check out?

        r = np.linspace(min_r, max_r, nr)
        t = np.linspace(min_t, max_t, nt)

        #--------------------------------------------
        # Define the cartesian grid
        #--------------------------------------------

        nx = nr
        ny = nr

        x = np.linspace(-max_r, max_r, nx)
        y = np.linspace(-max_r, max_r, ny)
        

        #--------------------------------------------
        # Initialise the 3D cartesian hdf5 file and its metadata
        # E3 field (q3d) corresponds to E3 field (cartesian)
        # E4 field (q3d) corresponds to E2 field (cartesian)
        #--------------------------------------------


        direction_shortcut = {
        '1' : 0,
        '3' : 1,
        '4' : 2
        }

        
        path_shortcut = file[f'mode_{mode_number}']['re'][direction_shortcut[dir]][0]

        
        self.prepare_hdf5_file(path_shortcut,target_directory ,min_z, max_z, max_r, nx, ny, nz) # <- check this and how it could be ammended
        
        #--------------------------------------------
        # Fill the 3D cartesian hdf5 file slice by slice
        #--------------------------------------------
        
        f = h5py.File(self.target_path, 'a')
        

        
        X, Y = np.meshgrid(x, y)
        
        R = np.sqrt(X**2 + Y**2)
        
        theta = np.arctan2(Y, X)

        
        

        print('before starting the loop')
        for i in range(nz - 1):
            
            # We build a transverse (r,t) slice
            #z = np.zeros([nr, nt])

            # We loop on all thetas
            #for j,tj in enumerate(t):
            
            interpolated_fields = self.compute_fields_for_slice(i, data, r, R, theta, mode_number)
            try:
                if dir=='1':
                    field = interpolated_fields['field_1']
                elif dir=='3':
                    field = interpolated_fields['field_2'] * np.cos(theta) - interpolated_fields['field_3'] * np.sin(theta)
                elif dir=='4':
                    field = interpolated_fields['field_2'] * np.sin(theta) + interpolated_fields['field_3'] * np.cos(theta)
            except IndexError:

                print(f'issue at {i}, {np.shape(interpolated_fields["field_1"])}')
                break    


            # We now interpolate Ex and Ey from the polar grid to the cartesian grid
            #simdata_cart = self.polar2cartesian(r, t, z, x, y, order=3)
            #print(np.shape(simdata_cart))
            f['/dataset'][:,:,i] = np.nan_to_num(field)

        
        f.close()

        return
