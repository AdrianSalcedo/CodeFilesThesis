# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gabrielsalcedo/Bocop-2.2.1-linux

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build

# Include any dependencies generated for this target.
include CMakeFiles/bocop.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bocop.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bocop.dir/flags.make

CMakeFiles/bocop.dir/core/main.cpp.o: CMakeFiles/bocop.dir/flags.make
CMakeFiles/bocop.dir/core/main.cpp.o: ../../core/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bocop.dir/core/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bocop.dir/core/main.cpp.o -c /home/gabrielsalcedo/Bocop-2.2.1-linux/core/main.cpp

CMakeFiles/bocop.dir/core/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bocop.dir/core/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gabrielsalcedo/Bocop-2.2.1-linux/core/main.cpp > CMakeFiles/bocop.dir/core/main.cpp.i

CMakeFiles/bocop.dir/core/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bocop.dir/core/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gabrielsalcedo/Bocop-2.2.1-linux/core/main.cpp -o CMakeFiles/bocop.dir/core/main.cpp.s

# Object files for target bocop
bocop_OBJECTS = \
"CMakeFiles/bocop.dir/core/main.cpp.o"

# External object files for target bocop
bocop_EXTERNAL_OBJECTS =

../bocop: CMakeFiles/bocop.dir/core/main.cpp.o
../bocop: CMakeFiles/bocop.dir/build.make
../bocop: lib/libCORE.a
../bocop: ../../ThirdParty/Ipopt-3.12.12/lib/libipopt.so
../bocop: ../../ThirdParty/Ipopt-3.12.12/lib/libcoinmumps.so
../bocop: ../../ThirdParty/Ipopt-3.12.12/lib/libcoinlapack.so
../bocop: ../../ThirdParty/Ipopt-3.12.12/lib/libcoinblas.so
../bocop: CMakeFiles/bocop.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bocop"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bocop.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bocop.dir/build: ../bocop

.PHONY : CMakeFiles/bocop.dir/build

CMakeFiles/bocop.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bocop.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bocop.dir/clean

CMakeFiles/bocop.dir/depend:
	cd /home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gabrielsalcedo/Bocop-2.2.1-linux /home/gabrielsalcedo/Bocop-2.2.1-linux /home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build /home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build /home/gabrielsalcedo/Bocop-2.2.1-linux/TomatoSystemControl/build/CMakeFiles/bocop.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bocop.dir/depend

