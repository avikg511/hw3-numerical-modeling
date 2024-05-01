# hw3-numerical-modeling
## Description
Here is a numerical model to solve the diffusion advection equation using Forward Euler and a mix of Forward Euler and Leap Frog/Centered Differencing. The math is largely explained, but follows conventions from the class.

The general use of the code is as follows. Currently, it's an XCode project so you should open it up with XCode. There are other branches with a CMakeSetup and this could theoretically be ported there without much work. Nevertheless, it's a lot simpler to debug on XCode as opposed to lldb (OSX gdb for CPP) so I switched to XCode. 

The way you work with the code is to build the code and run it. Then, using your terminal, you should run ./build.sh. Be sure to use `chmod +x ./build.sh` while in the directory with `./build.sh` to give it permissions to run.

Also, double check your scheme using Product >> Scheme >> Edit Scheme >> Options >> Set working directory. Otherwise, your code won't actually write anything to the CSV files it's supposed to. I tried building with GNUPlot from XCode, just using the system calls from CPP but it didn't work and I'm not sure why. Having a separate terminal instance is probably the simpler solution for now. 


## File Tree:
```
.
├── DiffusionModel.cpp
├── DiffusionModel.hpp
├── LICENSE
├── README.md
├── homework3
│   ├── FEData.csv
│   ├── FELFrogMixData.csv
│   ├── build.sh
│   └── main.cpp
└── homework3.xcodeproj
    ├── project.pbxproj
    ├── project.xcworkspace
    │   ├── contents.xcworkspacedata
    │   ├── xcshareddata
    │   │   ├── IDEWorkspaceChecks.plist
    │   │   └── swiftpm
    │   │       └── configuration
    │   └── xcuserdata
    │       └── aghos.xcuserdatad
    │           ├── UserInterfaceState.xcuserstate
    │           └── xcdebugger
    │               └── Expressions.xcexplist
    ├── xcshareddata
    │   └── xcschemes
    │       └── homework3.xcscheme
    └── xcuserdata
        └── aghos.xcuserdatad
            ├── xcdebugger
            │   └── Breakpoints_v2.xcbkptlist
            └── xcschemes
                └── xcschememanagement.plist

16 directories, 16 files

```

