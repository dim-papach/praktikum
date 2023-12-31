# maskfill

Fill in masked values in an image.

Method description is in van Dokkum et al. (2023) (PASP).

Examples in
 example_synthetic/
 example_m51/

## Installation 

The `maskfill` code can be installed from `PyPI` via 
```
pip install maskfill
```
or the latest version from github can be retrieved via 

```
git clone https://github.com/dokkum/maskfill.git
cd maskfill 
pip install . # or pip install -e .
```
Its dependencies are only `numpy`, `scipy`, `astropy` (for fits handling), and Python>=3.6 for type annotations. If you wish, you can create a minimal working environment for this code via `conda` or `mamba`, e.g.,  
```
mamba create -n maskfill python=3.10 numpy scipy astropy
activate maskfill
```

## Usage 

`Maskfill` can be used either from the command line via a command line interface (CLI), or imported in Python as a function. 

### CLI 
When you install `maskfill` it will create a `maskfill` executable callable from the shell. You can print the usage via 
```
maskfill -h
```

The simplest call is something like 
```
maskfill im.fits mask.fits out.fits
```

or
```
maskfill im mask out #file extension assumed to be .fits
```

in which we provide the input image, mask image, and name of the output file (if the `.fits` is omitted, `maskfill` will add it, though if your files have alternate extensions like `.fit` you should specify the full name). There are also several optional arguments and flags. 

- `-e X` or `--extension X`: if the image and mask are not in the 0th fits extension, specify it here 
- `-s X` or `--size X`: if you want a larger window kernel than the minimum 3x3, specify it here (faster, but less accurate results)
- `-o median` or `--operator median`: either 'median' or 'mean', defines how masked pixels are filled in based on their neighbors
- `-n` or `--nosmooth`: disable a final-step boxcar smoothing of the filled in mask pixels
- `-w` or `--writesteps`: write `_iter_N` fits files after each iteration of the algorithm (default is False)
- `-v` or `--verbose`: verbose output (shows the progress of iterations and number of remaining masked pixels).

The output is saved in the fits file with the provided output name. By default, after infilling, a smoothing step (using the same window, but a `mean` filter) is used to reduce sharp edges introduced by the iterative infilling. When enabled, the output fits file will contain the smoothed output image in the 0th extension, and the unsmoothed version post infilling in the 1st extension. If `nosmooth` is flagged, the 0th extension will contain the unsmoothed output. Information about which type of output is in which extension is added to the header. 

### Python Function 

In python, import the `maskfill` function via 
```python
from maskfill import maskfill
```
The function supports all the same features as the CLI. We can provide either file paths or arrays directly, e.g., 
```python
from astropy.io import fits 

im = fits.getdata('im.fits')
mask = fits.getdata('mask.fits')

result = maskfill(im,mask)
# or, equivalently, 
result = maskfill('im.fits','mask.fits')
```
When providing filepaths, the `ext` argument allows non 0th extensions to be queried. Like the CLI, all other arguments are optional, but unlike the CLI, the `output_file` argument is also optional --- `maskfill` returns the filled image(s) as arrays directly, but a fits file can be written if desired by providing a filename (same goes for the intermediate steps). 

Thus, a "maximally verbose" run might look like 
```python
result = maskfill('im.fits','mask.fits',writesteps=True,output_file='res.fits',verbose=True)
```
which would not only print to screen the iteration progress but write both intermediate and final fits files to disk. 

The `results` output by the function are a tuple containing either: 

1. `(output_smoothed,output)` if the smoothing step is enabled, or 
2. `(output, None)` if the smoothing step is disabled. 

Thus to access or plot the smoothed image from a default run, we would either run `results[0]`, or run with 

```python
results,_ = maskfill(im,mask)
```
