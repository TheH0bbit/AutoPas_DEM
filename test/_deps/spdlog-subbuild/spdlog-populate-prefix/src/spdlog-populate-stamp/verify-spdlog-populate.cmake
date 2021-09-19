# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip")
  message(FATAL_ERROR "File not found: /home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip")
endif()

if("MD5" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("7415a9768f3433bd93d78c1c87fd0576" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip'")

file("MD5" "/home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip" actual_value)

if(NOT "${actual_value}" STREQUAL "7415a9768f3433bd93d78c1c87fd0576")
  message(FATAL_ERROR "error: MD5 hash of
  /home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip
does not match expected value
  expected: '7415a9768f3433bd93d78c1c87fd0576'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
