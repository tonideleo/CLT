# **CLT - Classical Lamination Theory**

This is an object-oriented code written in Matlab as a tool for analyzing classical lamination theory problems.

So far it is composed by 3 main classes:

1)  **CLT_Material**
    This class will handle all the material aspect of the laminate.

2)  **CLT_Ply**
    All the operation at ply level (local stresses in X-Y or 1-2 directions) are done here.

3)  **CLT_Laminate**
    This class will manage and merge all the plies created with the two previous classes. 
    Multiple methods are available to help implementing books' equations and computational methods.

You can run the [test_input.m](https://github.com/tonideleo/CLT/blob/main/test_input.m) file to start. Here a quick description:

## MATERIAL

This is how you create an instance of the material class:
```
a   =   CLT_Material(E1         =   170e3,...
                     E2         =   10e3,...
                     G12        =   13e3,...
                     v12        =   0.3,...
                     Name       =   'IM7',...
                     Verbose    =   false);
```
Pretty self-explanatory, the Verbose option is not implemented yet.
You can add damage behavior (as of right now only Tsai-Wu is implemted) by method of the class:
```
a.addFailureTW(s1t       =   2.723e3,...
               s1c       =   1.2e3,...
               s2t       =   127,...
               s2c       =   200,...
               t12       =   95.8,...
               F12       =   0);
```
## LAMINATE
To create a laminate, you can create an instance of the laminate class in this way: 
```
b   =   CLT_Laminate(Stack          =   [0 90 90 0],...
                     Material       =   a,...
                     Thickness      =   0.15,...
                     Symmetric      =   false,...
                     RepeatLeft     =   1,...
                     RepeatRight    =   1,...
                     Verbose        =   true);
```                     
If plies have different materials (let's say a1, a2, a3, etc...) and different thicknesses, they can be explicitly written out:
```
b   =   CLT_Laminate(Stack          =   [0 45 -45 45],...
                     Material       =   [a1 a2 a3 a4],...
                     Thickness      =   [0.15 0.11 0.13 0.5],...
                     Symmetric      =   false,...
                     RepeatLeft     =   1,...
                     RepeatRight    =   1,...
                     Verbose        =   true);
```                 
The options Symmetric, RepeatLeft, and RepeatRight are used for normal conventions.
For example, if you have a laminate [0/45]2s2, instead of writing out
[0/45/0/45/45/0/45/0/0/45/0/45/45/0/45/0] you can the options set as:
```
b   =   CLT_Laminate(Stack          =   [0 45],...
                     Material       =   a,...
                     Thickness      =   0.15,...
                     Symmetric      =   true,...
                     RepeatLeft     =   2,...
                     RepeatRight    =   2,...
                     Verbose        =   true);
```               
You can request macroscopic effective behavior (E1,E2,G12,v12):
```
[Ex,Ey,Gxy,nu_xy] = b.calculateMacroBehavior() 
```
together with A,B,D,ABD,a,b,d,abd matrices:
```
A11 = b.A(1,1);
D66 = b.D(3,3); % or b.ABD(6,6)
% etc...
```
## Progressive Failure Analysis
If you use the following code:
```
b.ProgressiveFailureAnalysis(Nx    =   -100,...
                             Ny    =   0,...
                             Nxy   =   0,...
                             Mx    =   0,...
                             My    =   0,...
                             Mxy   =   0,...
                             Verbose   =   true);
```
you will obtain the following:
```
==================== Progressive Failure Analysis ====================
Ply 2 (Angle = 90 deg) FAILED at 68.56 % of given loads!
Ply 3 (Angle = 90 deg) FAILED at 68.56 % of given loads!
Ply 1 (Angle = 0 deg) FAILED at 81.51 % of given loads!
Ply 4 (Angle = 0 deg) FAILED at 81.51 % of given loads!
ALL plies failed!
Progressive Failure Analysis Successfully Terminated!
======================================================================
```
## Buckling Critical Load
Instead for buckling critical loading:
```
b.calculateBuckling(Lx      =   100,...
                    Ly      =   100,...
                    k       =   0,...
                    BC      =   'S4',...
                    DisplayPlot     =   true,...
                    DisplayTable    =   true);
```
with the following output:
```
    Critical Load    n    m
    _____________    _    _

       4.2457        1    1
       11.881        1    2
       15.376        2    1
       16.983        2    2
       25.201        1    3
       29.187        2    3
       30.844        3    2
       38.212        3    3
       55.357        3    1
```
Note: this method is not fully completed yet, as mostly the rest of the code.

Please feel free to reach out if you have any question or if you would like to have something implemented!

### Antonio Alessandro Deleo
### Department of Aeronautics & Astronautics
### University of Washington Seattle
### a(dot)deleo(at)uw.(edu)
