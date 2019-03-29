project(MD)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)


include_directories(
	${Boost_INCLUDE_DIR}
)

link_directories(
	${Boost_LIBRARYDIR}
)

set(COMPILER_FLAGS
	"-ggdb -Wall -L/usr/include -I/usr/include -std=c++0x -pipe -Wno-deprecated -pedantic"
)

add_executable(mod_dec
	Commons.h
	Token.h
	Token.cpp
	DyckWord.h
	DyckWord.cpp
	ModDecomp.h
	ModDecomp.cpp
)

set_target_properties(mod_dec PROPERTIES
	COMPILE_FLAGS ${COMPILER_FLAGS}
)

target_link_libraries(mod_dec ${Boost_LIBRARIES})


