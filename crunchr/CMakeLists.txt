cmake_minimum_required(VERSION 3.0)
project(crunchr)


find_package(htslib 1.4 COMPONENTS htslib)




add_library(smalt
  src/smalt/filio.c
  src/smalt/sequence.c
  src/smalt/elib.c
  src/smalt/array.c)

#include_directories(src)
#include_directories(src/smalt)

add_executable(desiguide src/desiguide.c)

#add_executable(crunchr crunchr.c)



target_link_libraries(desiguide
  smalt
  ${htslib_LIBRARIES})

#target_link_libraries(crunchr ${htslib_LIBRARIES})


install(TARGETS desiguide
	DESTINATION bin)
