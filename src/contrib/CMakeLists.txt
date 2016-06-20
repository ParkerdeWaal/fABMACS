enable_language(Fortran)
set(CONTRIB_PROGRAMS 
     #add here any programs you want to compile

)

foreach(PROG ${CONTRIB_PROGRAMS})
        add_executable(${PROG} ${PROG}.f ${NGMX_COMMON_SOURCE})
        target_link_libraries(${PROG} gfortran)
SET_TARGET_PROPERTIES(${PROG} PROPERTIES LINKER_LANGUAGE Fortran)
	set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG)

