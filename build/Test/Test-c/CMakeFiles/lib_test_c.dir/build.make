# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/etud/Bureau/HGS-CVRP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/etud/Bureau/HGS-CVRP/build

# Include any dependencies generated for this target.
include Test/Test-c/CMakeFiles/lib_test_c.dir/depend.make

# Include the progress variables for this target.
include Test/Test-c/CMakeFiles/lib_test_c.dir/progress.make

# Include the compile flags for this target's objects.
include Test/Test-c/CMakeFiles/lib_test_c.dir/flags.make

Test/Test-c/CMakeFiles/lib_test_c.dir/test.c.o: Test/Test-c/CMakeFiles/lib_test_c.dir/flags.make
Test/Test-c/CMakeFiles/lib_test_c.dir/test.c.o: ../Test/Test-c/test.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/etud/Bureau/HGS-CVRP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object Test/Test-c/CMakeFiles/lib_test_c.dir/test.c.o"
	cd /home/etud/Bureau/HGS-CVRP/build/Test/Test-c && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/lib_test_c.dir/test.c.o -c /home/etud/Bureau/HGS-CVRP/Test/Test-c/test.c

Test/Test-c/CMakeFiles/lib_test_c.dir/test.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lib_test_c.dir/test.c.i"
	cd /home/etud/Bureau/HGS-CVRP/build/Test/Test-c && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/etud/Bureau/HGS-CVRP/Test/Test-c/test.c > CMakeFiles/lib_test_c.dir/test.c.i

Test/Test-c/CMakeFiles/lib_test_c.dir/test.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lib_test_c.dir/test.c.s"
	cd /home/etud/Bureau/HGS-CVRP/build/Test/Test-c && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/etud/Bureau/HGS-CVRP/Test/Test-c/test.c -o CMakeFiles/lib_test_c.dir/test.c.s

# Object files for target lib_test_c
lib_test_c_OBJECTS = \
"CMakeFiles/lib_test_c.dir/test.c.o"

# External object files for target lib_test_c
lib_test_c_EXTERNAL_OBJECTS =

Test/Test-c/lib_test_c: Test/Test-c/CMakeFiles/lib_test_c.dir/test.c.o
Test/Test-c/lib_test_c: Test/Test-c/CMakeFiles/lib_test_c.dir/build.make
Test/Test-c/lib_test_c: libhgscvrp.so
Test/Test-c/lib_test_c: Test/Test-c/CMakeFiles/lib_test_c.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/etud/Bureau/HGS-CVRP/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable lib_test_c"
	cd /home/etud/Bureau/HGS-CVRP/build/Test/Test-c && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lib_test_c.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Test/Test-c/CMakeFiles/lib_test_c.dir/build: Test/Test-c/lib_test_c

.PHONY : Test/Test-c/CMakeFiles/lib_test_c.dir/build

Test/Test-c/CMakeFiles/lib_test_c.dir/clean:
	cd /home/etud/Bureau/HGS-CVRP/build/Test/Test-c && $(CMAKE_COMMAND) -P CMakeFiles/lib_test_c.dir/cmake_clean.cmake
.PHONY : Test/Test-c/CMakeFiles/lib_test_c.dir/clean

Test/Test-c/CMakeFiles/lib_test_c.dir/depend:
	cd /home/etud/Bureau/HGS-CVRP/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/etud/Bureau/HGS-CVRP /home/etud/Bureau/HGS-CVRP/Test/Test-c /home/etud/Bureau/HGS-CVRP/build /home/etud/Bureau/HGS-CVRP/build/Test/Test-c /home/etud/Bureau/HGS-CVRP/build/Test/Test-c/CMakeFiles/lib_test_c.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Test/Test-c/CMakeFiles/lib_test_c.dir/depend
