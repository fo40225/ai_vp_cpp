# CMakeList.txt: ai_vp_cpp 的 CMake 專案，在此包含來源及定義
# 專案專屬邏輯。
#

# 將來源新增到此專案的可執行檔。
add_executable (ai_vp_cpp "ai_vp_cpp.cpp" "ai_vp_cpp.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ai_vp_cpp PROPERTY CXX_STANDARD 20)
endif()

# TODO: 新增測試，並視需要安裝目標。
