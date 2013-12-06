#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

if (NOT DEFINED OUTPUT_DIR OR NOT DEFINED SOURCE_HTML_DIR)
    message(FATAL_ERROR "Required input parameter not set")
endif()

function(PRE_EXPORT_ACTIONS)
    file(MAKE_DIRECTORY ${OUTPUT_DIR})
    file(MAKE_DIRECTORY ${OUTPUT_DIR}/programs)
    file(COPY ${SOURCE_HTML_DIR}/header.html.in DESTINATION .)
    file(COPY ${SOURCE_HTML_DIR}/footer.html DESTINATION .)
endfunction()

function(POST_EXPORT_ACTIONS)
    # This is generated by gmx help -export html
    set(HEADER_FILE header.html)
    set(FOOTER_FILE ${SOURCE_HTML_DIR}/footer.html)
    file(READ ${HEADER_FILE} HEADER_TEXT)
    file(READ ${FOOTER_FILE} FOOTER_TEXT)
    set(_title_re "[Tt][Ii][Tt][Ll][Ee]")

    function(CREATE_HTML_FILE SOURCE_FILE ROOTPATH)
        file(RELATIVE_PATH _rel_path ${SOURCE_HTML_DIR} ${SOURCE_FILE})
        file(READ ${SOURCE_FILE} _content)
        string(REGEX REPLACE "^ *<${_title_re}>(.*)</${_title_re}>\n" "" _content "${_content}")
        set(TITLE "${CMAKE_MATCH_1}")
        string(CONFIGURE "${HEADER_TEXT}" _header @ONLY)
        set(_content "${_header}${_content}${FOOTER_TEXT}")
        file(WRITE ${OUTPUT_DIR}/${_rel_path} "${_content}")
    endfunction()

    create_html_file(${SOURCE_HTML_DIR}/online.html "")
    file(COPY ${SOURCE_HTML_DIR}/images DESTINATION ${OUTPUT_DIR})
    file(MAKE_DIRECTORY ${OUTPUT_DIR}/online)
    file(COPY ${SOURCE_HTML_DIR}/online/style.css DESTINATION ${OUTPUT_DIR}/online)
    file(GLOB _source_files ${SOURCE_HTML_DIR}/online/*.html)
    foreach(_file ${_source_files})
        create_html_file(${_file} "../")
    endforeach()
endfunction()

if (STEP STREQUAL "PRE")
    pre_export_actions()
elseif (STEP STREQUAL "POST")
    post_export_actions()
else()
    message(FATAL_ERROR "Unknown parameter STEP=${STEP}")
endif()
