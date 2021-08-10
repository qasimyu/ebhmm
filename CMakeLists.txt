cmake_minimum_required(VERSION 3.0)
project(EBHMM)

SET(CMAKE_BUILD_TYPE Release)

#add_definitions(-std=c++11)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory(lib)

include_directories (
	${EBHMM_SOURCE_DIR}/lib/bamtools/src
	${EBHMM_SOURCE_DIR}/lib/cnacaller
	${EBHMM_SOURCE_DIR}/lib/config
	${EBHMM_SOURCE_DIR}/lib/fastahack
	${EBHMM_SOURCE_DIR}/lib/genomedata
	${EBHMM_SOURCE_DIR}/lib/hmm
	${EBHMM_SOURCE_DIR}/lib/inputparser
	${EBHMM_SOURCE_DIR}/lib/kent/inc
	${EBHMM_SOURCE_DIR}/lib/mathfunc
	${EBHMM_SOURCE_DIR}/lib/matrix
	${EBHMM_SOURCE_DIR}/lib/model
	${EBHMM_SOURCE_DIR}/lib/mydefine
	${EBHMM_SOURCE_DIR}/lib/split
	${EBHMM_SOURCE_DIR}/lib/threadpool
)

add_executable(ebhmm src/ebhmm.cpp)

target_link_libraries(ebhmm mydefine)

SET_TARGET_PROPERTIES(ebhmm PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EBHMM_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET ebhmm PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET ebhmm PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(ebhmm "${CMAKE_THREAD_LIBS_INIT}")
endif()

INSTALL(TARGETS ebhmm
	RUNTIME DESTINATION bin
)
