
# Search for the maeparser library or clone the sources from GitHub
macro(find_maeparser maeparser_DIR MAEPARSER_SOURCE_TAG)

    find_package(maeparser QUIET)

    if(maeparser_FOUND)
        message(STATUS "Found compiled maeparser library at ${maeparser_DIR}")
    elseif(NOT maeparser_DIR STREQUAL "")
        message(FATAL_ERROR "*** Failed to find a compiled instance of maeparser under "
                    "${maeparser_DIR}.")
    else()
        find_package(Git QUIET)

        if(GIT_FOUND)
            message(STATUS "*** maeparser binary installation NOT found, getting/updating "
                    "maeparser sources from GitHub ...")

            if (EXISTS "${MAEPARSER_SRC_DIR}/maeparser/.git")
                execute_process(COMMAND ${GIT_EXECUTABLE} pull
                                WORKING_DIRECTORY ${MAEPARSER_SRC_DIR}/maeparser
                                RESULT_VARIABLE GIT_RESULT)
            else()
                execute_process(COMMAND ${GIT_EXECUTABLE} clone
                                https://github.com/schrodinger/maeparser.git
                                ${MAEPARSER_SRC_DIR}/maeparser WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                RESULT_VARIABLE GIT_RESULT)
            endif()

            if(NOT GIT_RESULT EQUAL "0")
                message(FATAL_ERROR "Failed to get maeparser from GitHub.")
            endif()

            execute_process(COMMAND ${GIT_EXECUTABLE} checkout ${MAEPARSER_SOURCE_TAG}
                            WORKING_DIRECTORY ${MAEPARSER_SRC_DIR}/maeparser
                            RESULT_VARIABLE GIT_RESULT)

            if(NOT GIT_RESULT EQUAL "0")
                message(FATAL_ERROR "Failed to check out tag ${MAEPARSER_SOURCE_TAG}.")
            endif()

        endif()

        if(NOT EXISTS "${MAEPARSER_SRC_DIR}/maeparser/CMakeLists.txt")
            message(FATAL_ERROR "Failed to find a valid instance of maeparser's "
                    "source code.")
        endif()

        set(maeparser_INCLUDE_DIRS ${MAEPARSER_SRC_DIR})
        add_subdirectory(${MAEPARSER_SRC_DIR}/maeparser)
    endif()

endmacro()