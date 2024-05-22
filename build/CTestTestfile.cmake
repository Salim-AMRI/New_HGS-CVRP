# CMake generated Testfile for 
# Source directory: /home/etud/Bureau/HGS-CVRP
# Build directory: /home/etud/Bureau/HGS-CVRP/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(bin_test_X-n101-k25 "/usr/bin/cmake" "-DINSTANCE=X-n101-k25" "-DCOST=27591" "-DROUND=1" "-P" "/home/etud/Bureau/HGS-CVRP/Test/TestExecutable.cmake")
set_tests_properties(bin_test_X-n101-k25 PROPERTIES  _BACKTRACE_TRIPLES "/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;32;add_test;/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;0;")
add_test(bin_test_X-n106-k14 "/usr/bin/cmake" "-DINSTANCE=X-n110-k13" "-DCOST=14971" "-DROUND=1" "-P" "/home/etud/Bureau/HGS-CVRP/Test/TestExecutable.cmake")
set_tests_properties(bin_test_X-n106-k14 PROPERTIES  _BACKTRACE_TRIPLES "/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;37;add_test;/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;0;")
add_test(bin_test_CMT6 "/usr/bin/cmake" "-DINSTANCE=CMT6" "-DCOST=555.43" "-DROUND=0" "-P" "/home/etud/Bureau/HGS-CVRP/Test/TestExecutable.cmake")
set_tests_properties(bin_test_CMT6 PROPERTIES  _BACKTRACE_TRIPLES "/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;44;add_test;/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;0;")
add_test(bin_test_CMT7 "/usr/bin/cmake" "-DINSTANCE=CMT7" "-DCOST=909.675" "-DROUND=0" "-P" "/home/etud/Bureau/HGS-CVRP/Test/TestExecutable.cmake")
set_tests_properties(bin_test_CMT7 PROPERTIES  _BACKTRACE_TRIPLES "/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;49;add_test;/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;0;")
add_test(lib_test_c "/home/etud/Bureau/HGS-CVRP/build/Test/Test-c/lib_test_c")
set_tests_properties(lib_test_c PROPERTIES  _BACKTRACE_TRIPLES "/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;61;add_test;/home/etud/Bureau/HGS-CVRP/CMakeLists.txt;0;")
subdirs("Test/Test-c")
