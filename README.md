# CS3353-Project4 - Solve A System of Linear Equations

## Overview
- The project mainly consists of two parts:
    - Compare Trivial (Gaussian Elimination) and Robust (Gauss-Seidel method): solve Ax = b where A, x, b are matrices
    - Introduce three new methods: Jacobi iterative method, Successive over-relaxation, and Thomas algorithm. Solve a more specific numeric problem: a second order differential equation. 
- The project does not require any command line argument if you want to run it from your terminal:
    <pre> python interface.py </pre>

## Menu instruction
- Press 1 to simply solve a system of linear equations(Gaussian Elimination & Gauss-Seidel); requires an input of path to the matrix file
- Press 2 to Run generate test (Gaussian Elimination vs Gauss-Seidel)
    - Under this option, you can input a directory that contains all the matrix files and a plot will be generated to show the performances.
- Press 3 to Run special numerical problem - second order ordinary differential equation
    - This option allows you to choose different methods among the five. Also you are allowed to change the data points to see different performances.
- Press E to exit
- instruction.txt is a sample execution (IDE: spyder through anaconda-navigator)

## Matrix file
- Matrix files are all .txt file with the size on the first line and a combined nx(n+1) matrix. File name represents the size n.
    <pre> Ex: 3.txt
    3
    4 1 2 4
    3 5 1 7
    1 1 3 3</pre>
    This represents the linear problem:
    <pre>
    [4 1 2] [x1]   [4]
    [3 5 1] [x2] = [7]
    [1 1 3] [x3]   [3]</pre>

## Some Info
- Project is written with python scripts. 
- Environment: &nbsp;conda 4.7.12  &nbsp; Python 3.7.4 &nbsp; IDE:  spyder==3.3.6 &nbsp;&nbsp;other dependencies are in requirements.txt<br>
- Video link:  &nbsp; https://www.dropbox.com/s/cdz2j3t4o7hnec7/Program4-demo.mp4?dl=0 <br>
- Video uses Jupyter Notebook for demo<br>
    - jupyter core     : 4.5.0<br>
    - jupyter-notebook : 6.0.1<br>
    - qtconsole        : 4.5.5<br>
    - ipython          : 7.8.0<br>
    - ipykernel        : 5.1.2<br>
    - jupyter client   : 5.3.3<br>
    - jupyter lab      : 1.1.4<br>
    - nbconvert        : 5.6.0<br>
    - ipywidgets       : 7.5.1<br>
    - nbformat         : 4.4.0<br>
    - traitlets        : 4.3.3<br>
    <br>
    - extra import for demo: <pre>import import_ipynb</pre>
        - installation: <pre> pip install import_ipynb </pre>
- /demo is a folder for jupyter notebooks, codes have extra demo functions
- /Graphs is a folder that contains test results with different data sets
- /matrices0 & /matrices1 & /small are all converging matrices that work for both Gaussian Elimination and Gauss-Seidel. It is also the test files.