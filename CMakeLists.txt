cmake_minimum_required(VERSION 3.14)


project(GLS CXX C)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 手动设置不同优化等级
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2") # 默认 -O2
# set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O1") # -O1
 set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0") # 无优化
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") # 强制添加 -O3

find_package(Boost 1.71.0 COMPONENTS program_options REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(
    GLS
    main.cpp
    src/readSeq.cpp
    src/fcgr.cpp
)

target_link_libraries(GLS PRIVATE
        ${Boost_LIBRARIES}
)