import numpy as np

def print_numpy_array_with_commas(array):
    # Check if the array is multi-dimensional
    if array.ndim > 1:
        # Handle multi-dimensional arrays
        subarrays = []
        for subarray in array:
            formatted_subarray = '[' + ', '.join(map(str, subarray)) + ']'
            subarrays.append(formatted_subarray)
        formatted_string = '[' + ', '.join(subarrays) + ']'
    else:
        # Handle one-dimensional array
        formatted_string = '[' + ', '.join(map(str, array)) + ']'
    
    print(formatted_string)

if __name__ == "__main__":
    python_datafiles = ["Tv_theta.npz", "Tv_theta_kxsource.npz", "Tv_theta_kxtarget.npz", "Tv_theta_kysource.npz"]

    for python_datafile in python_datafiles:
        transfer_python = np.load(python_datafile)
        transfer_python = transfer_python["T_v"]
        print("data for datafile " + python_datafile)
        print_numpy_array_with_commas(transfer_python)
        print("")
        print("")
        print("")