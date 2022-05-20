# R_IC user guide

Implementation of the Index of Connectivity in R language

Author details: Tommaso Baggio, Lorenzo Martini, Loris Torresani.
Script and data info: This script performs the Index of Connectivity (Cavalli et al., 2013).
Copyright statement: This script is the product of the work of Tommaso Baggio, Lorenzo Martini and Loris Torresani.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The release performs the Index of Connectivity (IC) starting from the version developed by Cavalli et al. (2013). We implemented IC in the R programming language to make it more accessible and versatile for final users. Furthermore, beyond the traditional implementation of computing IC, we developed two algorithms to compute IC recursively (batch process) and a more flexible version, allowing user-defined weighting factor and the possibility to perform an analysis of the outputs along a given longitudinal profile.

### Dependencies
In order to run the three scripts provided in the release, TauDEM has to be installed in order to successively call the executables. You can download and install the latest version of the complete windows installer from https://hydrology.usu.edu/taudem/taudem5/downloads.html To develop the three scripts, version number 5.3.7 was used, however, R_IC is expected to work also with earlier versions.
Since the script is written in R language, R software is required https://cran.r-project.org/bin/windows/base/ The adopted version is the 4.0.3. To compute operations, we further used built-in functions from different R packages (before installing and loading them remember to install Rtools 4 from https://cran.r-project.org/bin/windows/Rtools/rtools40.html). The necessary packages are raster (v 3.4.5), shapefiles (v 0.7), rgdal (v 1.5.19), cowplot (v 1.1.1) and ggplot2 (v 3.3.2). The reported versions are those used to develop the three scripts. However, we expect that earlier or updated versions of R and packages will be suitable as well.

### Usage 
Three scripts are available in the repository of R_IC. The scripts are structured in the same way: the first lines are reserved for the input files, followed by operations to check the input data. If the data are correct, the script continues with the core algorithm to compute IC and to save the outputs. The common inputs of the three scripts are:
•	Definition of the working directory: directory where input files are located and where output files will be saved. Remember to not use the symbol “/” as last character.
•	Names of the input files: the raster file representing the depitted DEM (the scripts do not perform any type of control on the input raster in terms of hydrological consistency) and the shapefile representing the target of the IC analysis. Note that these files need to be in the same reference coordinate system. The formats of the input files are those accepted by the functions raster and readOGR of the package raster and rgdal, respectively. In the script, batch_R_IC remember to not specify the extension of the shapefile (as reported in the training data).
•	Parameters representing (i) the size of the moving window to derive the roughness map and then the weighting factor expressed as cell number (only odd numbers greater than 1 are acceptable) and (ii) a flag to save intermediate outputs (flag = 1) or to keep only essential outputs as IC, Ddown, Dup and roughness (flag=0).
