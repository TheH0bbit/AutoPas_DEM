# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.16.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(spdlog-populate NONE)

include(ExternalProject)
ExternalProject_Add(spdlog-populate
                     "UPDATE_DISCONNECTED" "False" "URL" "/home/raphael/Bachelor/repo/AutoPas/libs/spdlog-1.x.zip" "URL_HASH" "MD5=7415a9768f3433bd93d78c1c87fd0576"
                    SOURCE_DIR          "/home/raphael/Bachelor/repo/AutoPas/test/_deps/spdlog-src"
                    BINARY_DIR          "/home/raphael/Bachelor/repo/AutoPas/test/_deps/spdlog-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)
