cmake_minimum_required (VERSION 3.0)


# options for MAGI
option (USE_HDF5 "On to build support for HDF5." ON)
option (USE_MPI "On to build support for MPI." OFF)
option (USE_SFMT "On to build support for SFMT." ON)
option (USE_SFMT19937 "On to use SFMT of period 2^19937 - 1." ON)
option (USE_SFMT-JUMP "On to build support for SFMT-jump." ON)
option (USE_TIPSY_FORMAT "On to output particle data in TIPSY format." OFF)
option (USE_GOTHIC "On to use the output for GOTHIC." OFF)
option (USE_GOTHIC_BLOCK "On to enable the block time step in GOTHIC." ON)
option (QUICK_CHECK "On to output partial particle data in ASCII format." OFF)


# global setting(s)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# if HDF5 is used, then check whether the parallel HDF5 or the serial HDF5
if (USE_HDF5)
  find_program (PHDF5_CONFIG "h5pcc")
  if (EXISTS ${PHDF5_CONFIG})
    set (USE_PHDF5 ON)
  else (EXISTS ${PHDF5_CONFIG})
    set (USE_PHDF5 OFF)
  endif (EXISTS ${PHDF5_CONFIG})
endif (USE_HDF5)


# if the parallel HDF5 is set, then MPI compiler is required.
if (USE_PHDF5)
  set (USE_MPI ON)
  # use mpicc if available
  find_program (MPICC "mpicc")
  if (MPICC)
    set (CMAKE_C_COMPILER "mpicc")
  endif (MPICC)
endif (USE_PHDF5)


# when USE_MPI is off, then use icc if available
if (NOT USE_MPI)
  find_program (INTELCC "icc")
  if (INTELCC)
    set (CMAKE_C_COMPILER "icc")
  endif (INTELCC)
endif (NOT USE_MPI)


# compile target
project (magi C)
set (EXECUTABLE_OUTPUT_PATH bin)
set (CMAKE_VERBOSE_MAKEFILE TRUE)
add_definitions (-DNDEBUG)
file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/dat)
file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc)
file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/log)
file (COPY samples/cfg DESTINATION ${CMAKE_BINARY_DIR})
file (COPY samples/sh DESTINATION ${CMAKE_BINARY_DIR})
file (COPY samples/py DESTINATION ${CMAKE_BINARY_DIR})


# source files
include_directories (src)
set (CCFILES src/init/magi.c)
set (CCFILES ${CCFILES} src/init/abel.c src/init/eddington.c src/init/king.c src/init/profile.c src/init/spline.c src/init/table.c)
set (CCFILES ${CCFILES} src/init/blas.c src/init/diskDF.c src/init/potdens.c)
set (CCFILES ${CCFILES} src/misc/allocate.c)
set (CCFILES ${CCFILES} src/file/io.c)


# library files
include_directories (inc)
set (CCLIBS inc/constants.c inc/rand.c inc/rotate.c inc/timer.c inc/myutil.c)


# if mpicc is required and not available, then search MPI
if (USE_MPI)
  set (CCLIBS ${CCLIBS} inc/mpilib.c)
  if (MPICC-NOTFOUND)
    find_package (MPI REQUIRED)
    include_directories ("${MPI_INCLUDE_PATH}")
    target_link_libraries (${PROJECT_NAME} "${MPI_LIBRARIES}")
  endif (MPICC-NOTFOUND)
else (USE_MPI)
  add_definitions (-DDISABLE_MPI)
endif (USE_MPI)


# flags for warning and performance optimization
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
if ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcheck -Wdeprecated -Wformat -Winline -Wsign-compare")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-aliasing -Wstrict-prototypes -Wtrigraphs -Wwrite-strings")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wuninitialized -Wunknown-pragmas -Wunused-function -Wunused-variable")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ip -xHost")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -opt-report 3 -vec-report3")
else ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra -Wunused -Wcast-qual -Wcast-align -Wmissing-declarations")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math -funroll-loops -march=native")
endif ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")


# flag to activate OpenMP
find_package (OpenMP REQUIRED)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
if (NOT DEFINED OMP_THREADS)
  cmake_host_system_information (RESULT PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)
  set (OMP_THREADS ${PHYSICAL_CORES})
endif (NOT DEFINED OMP_THREADS)
add_definitions ("-DCPUS_PER_PROCESS=(${OMP_THREADS})")


# include SFMT (optional)
if (USE_SFMT)
  if (NOT DEFINED SFMT_DIR)
    set (SFMT_DIR "${CMAKE_SOURCE_DIR}/SFMT")
  endif (NOT DEFINED SFMT_DIR)

  set (SFMT_SRC ${SFMT_DIR}/SFMT.c)
  if (EXISTS ${SFMT_SRC})
    add_definitions (-DUSE_SFMT -DHAVE_SSE2)
    include_directories (${SFMT_DIR})
    set (SFMT_CCFILES ${SFMT_SRC})
    if (USE_SFMT19937)
      add_definitions (-DSFMT_MEXP=19937)
    else (USE_SFMT19937)
      add_definitions (-DSFMT_MEXP=216091)
    endif (USE_SFMT19937)
  else (EXISTS ${SFMT_SRC})
    set (USE_SFMT OFF)
    set (USE_SFMT19937 OFF)
    set (USE_SFMT-JUMP OFF)
  endif (EXISTS ${SFMT_SRC})

  # include SFMT-jump (optional)
  if (USE_SFMT-JUMP)
    if (NOT DEFINED SFMT-JUMP_DIR)
      set (SFMT-JUMP_DIR "${CMAKE_SOURCE_DIR}/SFMT-Jump")
    endif (NOT DEFINED SFMT-JUMP_DIR)

    set (SFMT-JUMP_SRC ${SFMT-JUMP_DIR}/SFMT-jump.c)
    if (EXISTS ${SFMT-JUMP_SRC})
      add_definitions (-DUSE_SFMTJUMP)
      include_directories(${SFMT-JUMP_DIR})
      set (SFMT_CCFILES ${SFMT_CCFILES} ${SFMT-JUMP_SRC} inc/sfmtjump_polynomial.h)
    else (EXISTS ${SFMT-JUMP_SRC})
      set (USE_SFMT-JUMP OFF)
    endif (EXISTS ${SFMT-JUMP_SRC})
  endif (USE_SFMT-JUMP)

endif (USE_SFMT)


# include HDF5 (optional)
if (USE_HDF5)
  set (CCLIBS ${CCLIBS} inc/hdf5lib.c)
  find_package (HDF5 REQUIRED)
  add_definitions (-DUSE_HDF5_FORMAT)
  include_directories ("${HDF5_INCLUDE_DIRS}")
else (USE_HDF5)
  set (USE_TIPSY_FORMAT ON)
endif (USE_HDF5)


# include GSL (required)
find_package (GSL REQUIRED)
include_directories ("${GSL_INCLUDE_DIRS}")


# connect with GOTHIC (required)
if (USE_GOTHIC)
  set (CCFILES ${CCFILES} src/misc/brent.c src/misc/tune.c)
  add_definitions (-DLOCALIZE_I_PARTICLES -DUSE_BRENT_METHOD)
  if (USE_HDF5)
    add_definitions (-DMONITOR_ENERGY_ERROR)
  endif (USE_HDF5)
  if (USE_GOTHIC_BLOCK)
    add_definitions (-DBLOCK_TIME_STEP)
  endif (USE_GOTHIC_BLOCK)
else (USE_GOTHIC)
  set (USE_GOTHIC_BLOCK OFF)
  add_definitions (-DRUN_WITHOUT_GOTHIC)
endif (USE_GOTHIC)


# use TIPSY format (optional)
if (USE_TIPSY_FORMAT)
  add_definitions (-DWRITE_IN_TIPSY_FORMAT)
endif (USE_TIPSY_FORMAT)


# print out ASCII data for quick check (optional)
if (QUICK_CHECK)
  add_definitions (-DPRINT_OUT_ASCII_DATA_FOR_QUICK_CHECK)
endif (QUICK_CHECK)


add_executable (${PROJECT_NAME} ${CCFILES} ${CCLIBS} ${SFMT_CCFILES})


# link HDF5
if (USE_HDF5)
  if (INTELCC)
    list(REMOVE_ITEM HDF5_LIBRARIES /usr/lib64/libm.so;/usr/lib/libm.so;/usr/lib32/libm.so)
  endif (INTELCC)
  target_link_libraries (${PROJECT_NAME} "${HDF5_LIBRARIES}")
endif (USE_HDF5)


# link GSL
if (INTELCC)
  # remove -lm
  string (REPLACE " -lm" "" GSL_LIBRARIES ${GSL_LIBRARIES})
endif (INTELCC)
target_link_libraries (${PROJECT_NAME} "${GSL_LIBRARIES}")


add_executable (sample src/init/sample.c inc/macro.h)
if (NOT INTELCC)
  target_link_libraries (sample "-lm")
endif (NOT INTELCC)


# compile another target
add_executable (read src/init/read.c src/file/io.c ${CCLIBS})
if (USE_HDF5)
  target_link_libraries (read "${HDF5_LIBRARIES}")
endif (USE_HDF5)
target_link_libraries (read "${GSL_LIBRARIES}")


# print out the compile options
message (STATUS "Generate ${PROJECT_NAME} ...")
message (STATUS "OMP_THREADS = ${OMP_THREADS}")
message (STATUS "USE_HDF5 = ${USE_HDF5}")
if (USE_HDF5)
  message (STATUS "USE_PHDF5 = ${USE_PHDF5}")
endif (USE_HDF5)
message (STATUS "USE_MPI = ${USE_MPI}")
message (STATUS "USE_TIPSY_FORMAT = ${USE_TIPSY_FORMAT}")
message (STATUS "USE_SFMT = ${USE_SFMT}")
if (USE_SFMT)
  message (STATUS "SFMT_DIR = ${SFMT_DIR}")
  message (STATUS "USE_SFMT19937 = ${USE_SFMT19937}")
  message (STATUS "USE_SFMT-JUMP = ${USE_SFMT-JUMP}")
  if (USE_SFMT-JUMP)
    message (STATUS "SFMT-JUMP_DIR = ${SFMT-JUMP_DIR}")
  endif (USE_SFMT-JUMP)
endif (USE_SFMT)
message (STATUS "USE_GOTHIC = ${USE_GOTHIC}")
if (USE_GOTHIC)
  message (STATUS "USE_GOTHIC_BLOCK = ${USE_GOTHIC_BLOCK}")
endif (USE_GOTHIC)
message (STATUS "QUICK_CHECK = ${QUICK_CHECK}")
