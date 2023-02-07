# Project-specific settings
PROJECT := evo_comp_ecology
EMP_DIR := ../../EmpiricalSafe/include/emp
PYBIND_DIR := $(shell python3 -m pybind11 --includes)

INCLUDE := $(shell python3 -m pybind11 --includes)
FLAG := -Wall -shared -std=c++17 -fPIC -fvisibility=hidden
SUFFIX := $(shell python3-config --extension-suffix)

# Flags to use regardless of compiler
CFLAGS_all := -Wall -shared -std=c++1z -fPIC -Wno-unused-function -I$(EMP_DIR)/ $(PYBIND_DIR)/ -I/usr/include/python3.8/ -o example$(shell python3-config --extension-suffix)

# Native compiler information
CXX := clang++
CXX_nat := $(CXX)
CFLAGS_nat := -O3 -DNDEBUG $(CFLAGS_all)
CFLAGS_nat_debug := -fprofile-arcs -ftest-coverage -g $(CFLAGS_all)

# Emscripten compiler information
# CXX_web := emcc
# OFLAGS_web_all := -s TOTAL_MEMORY=67108864 --js-library $(EMP_DIR)/web/library_emp.js --js-library $(EMP_DIR)/web/d3/library_d3.js -s EXPORTED_FUNCTIONS="['_main', '_empCppCallback']" -s DISABLE_EXCEPTION_CATCHING=1 -s NO_EXIT_RUNTIME=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['cwrap', 'stringToUTF8', 'writeStringToMemory']"#--embed-file configs
# OFLAGS_web := -Oz -DNDEBUG
# OFLAGS_web_debug := -Oz -pedantic -Wno-dollar-in-identifier-extension

# CFLAGS_web := $(CFLAGS_all) $(OFLAGS_web) $(OFLAGS_web_all)
# CFLAGS_web_debug := $(CFLAGS_all) $(OFLAGS_web_debug) $(OFLAGS_web_all)

default: test

test: tests/test_interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat_debug) tests/test_interaction_networks.cc 
#	./test


clean:
	rm -f *~ source/*.o test *.gcda *.gcno

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
