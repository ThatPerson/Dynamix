cmake_minimum_required(VERSION 3.17)
project(dynamix C)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})



set(CMAKE_C_STANDARD 99)

include_directories(./src)
include_directories(./src/models)

find_package(OpenMP)

add_executable(dynamix
        src/main.c
        src/runners.c
        src/errors.c
        src/read_data.c
        src/global_gaf.c
        src/datatypes.c
        src/verification.c
        src/chisq.c
        src/crosen.c
        src/anneal.c
        src/models/model.c)

add_executable(prate
        src/predict_rate.c
        src/runners.c
        src/errors.c
        src/read_data.c
        src/global_gaf.c
        src/datatypes.c
        src/verification.c
        src/chisq.c
        src/crosen.c
        src/anneal.c
        src/models/model.c)

add_executable(tests
        src/test/tests.c
        src/runners.c
        src/errors.c
        src/read_data.c
        src/global_gaf.c
        src/datatypes.c
        src/verification.c
        src/chisq.c
        src/crosen.c
        src/anneal.c
        src/models/model.c
        src/test/old_model/smf.c
        src/test/old_model/emf.c
        src/test/old_model/aimf.c
        src/test/old_model/egaf.c
        src/test/old_model/gaf.c )

add_executable(correl
        src/correlation.c
        src/runners.c
        src/errors.c
        src/read_data.c
        src/global_gaf.c
        src/datatypes.c
        src/verification.c
        src/chisq.c
        src/crosen.c
        src/anneal.c
        src/models/model.c)

target_link_libraries(dynamix ${MPI_LIBRARIES})
target_link_libraries(dynamix m)
target_link_libraries(prate ${MPI_LIBRARIES})
target_link_libraries(prate m)
target_link_libraries(tests ${MPI_LIBRARIES})
target_link_libraries(tests m)
target_link_libraries(tests ${MPI_LIBRARIES})

target_link_libraries(correl m)
target_link_libraries(correl ${MPI_LIBRARIES})
target_link_libraries(tests cmocka)

set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -fopenmp ")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Wextra -fopenmp -fprofile-arcs -ftest-coverage")

if("${CMAKE_C_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang"
        OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang")
    message("Building debug with llvm Code Coverage")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fprofile-instr-generate -fcoverage-mapping")
elseif(CMAKE_COMPILER_IS_GNUC)
    message("Building debug with lcov Code Coverage Tools")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
endif()

if(MPI_COMPILE_FLAGS)
    set_target_properties(dynamix PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(tests PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(dynamix PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(tests PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")

endif()
