# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.4

# Default target executed when no arguments are given to make.
default_target: all

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ted/work/boone/ccqe/xs_ted/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ted/work/boone/ccqe/xs_ted

# Include the progress variables for this target.
include CMakeFiles/progress.make

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/ted/work/boone/ccqe/xs_ted/CMakeFiles $(CMAKE_ALL_PROGRESS)
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/ted/work/boone/ccqe/xs_ted/CMakeFiles 0

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean

# The main clean target
clean/fast: clean

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1

#=============================================================================
# Target rules for targets named foo

# Build rule for target.
foo: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 foo

# fast build rule for target.
foo/fast:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/build

# target to build an object file
src/TT_drawer.o:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_drawer.o

# target to preprocess a source file
src/TT_drawer.i:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_drawer.i

# target to generate assembly for a file
src/TT_drawer.s:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_drawer.s

# target to build an object file
src/TT_event.o:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_event.o

# target to preprocess a source file
src/TT_event.i:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_event.i

# target to generate assembly for a file
src/TT_event.s:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_event.s

# target to build an object file
src/TT_generator.o:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_generator.o

# target to preprocess a source file
src/TT_generator.i:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_generator.i

# target to generate assembly for a file
src/TT_generator.s:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_generator.s

# target to build an object file
src/TT_nucleus.o:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_nucleus.o

# target to preprocess a source file
src/TT_nucleus.i:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_nucleus.i

# target to generate assembly for a file
src/TT_nucleus.s:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_nucleus.s

# target to build an object file
src/TT_params.o:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_params.o

# target to preprocess a source file
src/TT_params.i:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_params.i

# target to generate assembly for a file
src/TT_params.s:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/TT_params.s

# target to build an object file
src/foo.o:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/foo.o

# target to preprocess a source file
src/foo.i:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/foo.i

# target to generate assembly for a file
src/foo.s:
	$(MAKE) -f CMakeFiles/foo.dir/build.make CMakeFiles/foo.dir/src/foo.s

# Help Target
help::
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... foo"
	@echo "... rebuild_cache"
	@echo "... src/TT_drawer.o"
	@echo "... src/TT_drawer.i"
	@echo "... src/TT_drawer.s"
	@echo "... src/TT_event.o"
	@echo "... src/TT_event.i"
	@echo "... src/TT_event.s"
	@echo "... src/TT_generator.o"
	@echo "... src/TT_generator.i"
	@echo "... src/TT_generator.s"
	@echo "... src/TT_nucleus.o"
	@echo "... src/TT_nucleus.i"
	@echo "... src/TT_nucleus.s"
	@echo "... src/TT_params.o"
	@echo "... src/TT_params.i"
	@echo "... src/TT_params.s"
	@echo "... src/foo.o"
	@echo "... src/foo.i"
	@echo "... src/foo.s"



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0

