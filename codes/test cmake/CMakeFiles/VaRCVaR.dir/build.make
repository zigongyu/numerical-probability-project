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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake"

# Include any dependencies generated for this target.
include CMakeFiles/VaRCVaR.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/VaRCVaR.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/VaRCVaR.dir/flags.make

CMakeFiles/VaRCVaR.dir/src/main.cpp.o: CMakeFiles/VaRCVaR.dir/flags.make
CMakeFiles/VaRCVaR.dir/src/main.cpp.o: src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/VaRCVaR.dir/src/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VaRCVaR.dir/src/main.cpp.o -c "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake/src/main.cpp"

CMakeFiles/VaRCVaR.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VaRCVaR.dir/src/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake/src/main.cpp" > CMakeFiles/VaRCVaR.dir/src/main.cpp.i

CMakeFiles/VaRCVaR.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VaRCVaR.dir/src/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake/src/main.cpp" -o CMakeFiles/VaRCVaR.dir/src/main.cpp.s

# Object files for target VaRCVaR
VaRCVaR_OBJECTS = \
"CMakeFiles/VaRCVaR.dir/src/main.cpp.o"

# External object files for target VaRCVaR
VaRCVaR_EXTERNAL_OBJECTS =

VaRCVaR: CMakeFiles/VaRCVaR.dir/src/main.cpp.o
VaRCVaR: CMakeFiles/VaRCVaR.dir/build.make
VaRCVaR: lib/libVaRCVaRlib.a
VaRCVaR: CMakeFiles/VaRCVaR.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable VaRCVaR"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/VaRCVaR.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/VaRCVaR.dir/build: VaRCVaR

.PHONY : CMakeFiles/VaRCVaR.dir/build

CMakeFiles/VaRCVaR.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/VaRCVaR.dir/cmake_clean.cmake
.PHONY : CMakeFiles/VaRCVaR.dir/clean

CMakeFiles/VaRCVaR.dir/depend:
	cd "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake" "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake" "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake" "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake" "/Users/zigong/PetF/S2/projet info/numerical-probability-project/codes/test cmake/CMakeFiles/VaRCVaR.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/VaRCVaR.dir/depend
