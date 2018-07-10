project(tables)
pkg_check_modules(cfitsio cfitsio>=3.2.8 REQUIRED)
add_library(${PROJECT_NAME}
	addCol.c
	c_imt.c
	c_irafinit.c
	c_phopar.c
	c_tbagt.c
	c_tbapt.c
	c_tbcdef1.c
	c_tbcfnd1.c
	c_tbciga.c
	c_tbcigi.c
	c_tbcigt.c
	c_tbcinf.c
	c_tbcnum.c
	c_tbegt.c
	c_tbept.c
	c_tbfpri.c
	c_tbfres.c
	c_tbhad.c
	c_tbhcal.c
	c_tbhgcm.c
	c_tbhgnp.c
	c_tbhgt.c
	c_tbhpcm.c
	c_tbhptt.c
	c_tbparse.c
	c_tbpsta.c
	c_tbrcsc.c
	c_tbrgtr.c
	c_tbrudf.c
	c_tbtacc.c
	c_tbtclo.c
	c_tbtcre.c
	c_tbtnam.c
	c_tbtopn.c
	c_vfn2osfn.c
	cerror.c
	ctables.c
	initCol.c
	tbCopyTmpl.c
	tbHeader.c
	tbSaveInfo.c
	tbl_util.c
)
include_directories(${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME}
	${cfitsio_LIBRARIES}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${cfitsio_INCLUDE_DIRS}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_compile_options(${PROJECT_NAME}
	PUBLIC ${cfitsio_CFLAGS_OTHER}
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)

