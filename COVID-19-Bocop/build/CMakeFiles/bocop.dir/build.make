# CMAKE generated file: DO NOT EDIT!
# Generated by "MSYS Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "/C/Program Files (x86)/CMake/bin/cmake.exe"

# The command to remove a file.
RM = "/C/Program Files (x86)/CMake/bin/cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /D/Bocop_Determinista/Bocop-2.2.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/build

# Include any dependencies generated for this target.
include CMakeFiles/bocop.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bocop.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bocop.dir/flags.make

CMakeFiles/bocop.dir/core/main.cpp.obj: CMakeFiles/bocop.dir/flags.make
CMakeFiles/bocop.dir/core/main.cpp.obj: CMakeFiles/bocop.dir/includes_CXX.rsp
CMakeFiles/bocop.dir/core/main.cpp.obj: D:/Bocop_Determinista/Bocop-2.2.1/core/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bocop.dir/core/main.cpp.obj"
	/D/BocopHJB-1.1.0/MinGW/bin/g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bocop.dir/core/main.cpp.obj -c /D/Bocop_Determinista/Bocop-2.2.1/core/main.cpp

CMakeFiles/bocop.dir/core/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bocop.dir/core/main.cpp.i"
	/D/BocopHJB-1.1.0/MinGW/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /D/Bocop_Determinista/Bocop-2.2.1/core/main.cpp > CMakeFiles/bocop.dir/core/main.cpp.i

CMakeFiles/bocop.dir/core/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bocop.dir/core/main.cpp.s"
	/D/BocopHJB-1.1.0/MinGW/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /D/Bocop_Determinista/Bocop-2.2.1/core/main.cpp -o CMakeFiles/bocop.dir/core/main.cpp.s

CMakeFiles/bocop.dir/core/main.cpp.obj.requires:

.PHONY : CMakeFiles/bocop.dir/core/main.cpp.obj.requires

CMakeFiles/bocop.dir/core/main.cpp.obj.provides: CMakeFiles/bocop.dir/core/main.cpp.obj.requires
	$(MAKE) -f CMakeFiles/bocop.dir/build.make CMakeFiles/bocop.dir/core/main.cpp.obj.provides.build
.PHONY : CMakeFiles/bocop.dir/core/main.cpp.obj.provides

CMakeFiles/bocop.dir/core/main.cpp.obj.provides.build: CMakeFiles/bocop.dir/core/main.cpp.obj


# Object files for target bocop
bocop_OBJECTS = \
"CMakeFiles/bocop.dir/core/main.cpp.obj"

# External object files for target bocop
bocop_EXTERNAL_OBJECTS =

C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: CMakeFiles/bocop.dir/core/main.cpp.obj
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: CMakeFiles/bocop.dir/build.make
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: lib/libCORE.a
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: D:/Bocop_Determinista/Bocop-2.2.1/ThirdParty/Ipopt-3.12.12/lib/libipopt.a
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: D:/Bocop_Determinista/Bocop-2.2.1/ThirdParty/Ipopt-3.12.12/lib/libcoinmumps.a
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: D:/Bocop_Determinista/Bocop-2.2.1/ThirdParty/Ipopt-3.12.12/lib/libcoinlapack.a
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: D:/Bocop_Determinista/Bocop-2.2.1/ThirdParty/Ipopt-3.12.12/lib/libcoinblas.a
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: CMakeFiles/bocop.dir/linklibs.rsp
C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe: CMakeFiles/bocop.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe"
	"/C/Program Files (x86)/CMake/bin/cmake.exe" -E remove -f CMakeFiles/bocop.dir/objects.a
	/D/BocopHJB-1.1.0/MinGW/bin/ar.exe cr CMakeFiles/bocop.dir/objects.a @CMakeFiles/bocop.dir/objects1.rsp
	/D/BocopHJB-1.1.0/MinGW/bin/g++.exe -std=gnu++11 -O -DNDEBUG -w   -Wl,--whole-archive CMakeFiles/bocop.dir/objects.a -Wl,--no-whole-archive  -o /C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe -Wl,--out-implib,/C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/libbocop.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/bocop.dir/linklibs.rsp

# Rule to build all files generated by this target.
CMakeFiles/bocop.dir/build: C:/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/bocop.exe

.PHONY : CMakeFiles/bocop.dir/build

CMakeFiles/bocop.dir/requires: CMakeFiles/bocop.dir/core/main.cpp.obj.requires

.PHONY : CMakeFiles/bocop.dir/requires

CMakeFiles/bocop.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bocop.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bocop.dir/clean

CMakeFiles/bocop.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MSYS Makefiles" /D/Bocop_Determinista/Bocop-2.2.1 /D/Bocop_Determinista/Bocop-2.2.1 /C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/build /C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/build /C/Users/Usuario1/Desktop/Carpeta/COVID-19-Bocop/build/CMakeFiles/bocop.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bocop.dir/depend

