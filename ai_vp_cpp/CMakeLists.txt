# CMakeList.txt: ai_vp_cpp 的 CMake 專案，在此包含來源及定義
# 專案專屬邏輯。
#

# 將來源新增到此專案的可執行檔。
add_executable (ai_vp_cpp "ai_vp_cpp.cpp" "ai_vp_cpp.h")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    find_package(simdjson REQUIRED)
    target_link_libraries(ai_vp_cpp PRIVATE simdjson::simdjson)
else()
    include(FetchContent)

    FetchContent_Declare(
        simdjson
        GIT_REPOSITORY https://github.com/simdjson/simdjson.git
        GIT_TAG  tags/v3.6.3
        GIT_SHALLOW TRUE)
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_MakeAvailable(simdjson)

    target_link_libraries(ai_vp_cpp PRIVATE simdjson)
endif()

find_package(Boost REQUIRED)
target_include_directories(ai_vp_cpp PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ai_vp_cpp PRIVATE ${Boost_LIBRARIES})

find_package(nlohmann_json REQUIRED)
target_link_libraries(ai_vp_cpp PRIVATE nlohmann_json::nlohmann_json)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ai_vp_cpp PROPERTY CXX_STANDARD 20)
endif()

# TODO: 新增測試，並視需要安裝目標。
