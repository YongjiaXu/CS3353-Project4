# CS3353-Project4 - Solve A System of Linear Equations

## Overview
- The project mainly consists of two parts:
    - Compare Trivial (Gaussian Elimination) and Robust (Gauss-Seidel method): solve Ax = b where A, x, b are matrices
    - Introduce three new methods: Jacobi iterative method, Successive over-relaxation, and Thomas algorithm. Solve a more specific numeric problem: a second order differential equation. 
- The project has two modes:
    - Simply solve Ax = b using Gaussian Elimination and Gauss-Seidel method. It takes a command line argument of one single file that contains the matrix. Solutions and runtimes will be printed to the console.
    <pre> Ex: python interface.py 10_0.txt</pre>
    - Interactive menu mode where you can make more customized changes and compare the performance
    <pre> Ex: python interface.py </pre>

## Menu instruction
- Press 1 to Run generate test (Gaussian Elimination vs Gauss-Seidel)
    - Under this option, you can input a directory that contains all the matrix file and a plot will be generated to show the performances.
- Press 2 to Run special numerical problem - second order ordinary differential equation
    - This option allows you to choose different methods among the five. Also you are allowed to change the data points to see different performances.

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
    [1 1 3] [x3]   [3]
    </pre>


- Video demo is run on JupyterNote book; .ipynb files are just for demo purpose
- Environment: conda 4.7.12 Python 3.7.4 spyder==3.3.6
- Jupyter env
  jupyter core     : 4.5.0
  jupyter-notebook : 6.0.1
  qtconsole        : 4.5.5
  ipython          : 7.8.0
  ipykernel        : 5.1.2
  jupyter client   : 5.3.3
  jupyter lab      : 1.1.4
  nbconvert        : 5.6.0
  ipywidgets       : 7.5.1
  nbformat         : 4.4.0
  traitlets        : 4.3.3
