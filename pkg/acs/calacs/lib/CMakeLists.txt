project(acs)
add_library(${PROJECT_NAME} SHARED
	acshist.c
	acsinfo.c
	acssect.c
	addk2d.c
	bin2d.c
	bincoords.c
	binupdate.c
	comparenum.c
	defswitch.c
	detchip.c
	detnsegn.c
	div1d.c
	dodqi.c
	donoise.c
	dostat.c
	err.c
	fileexists.c
	findbin.c
	findroot.c
	fromlt.c
	getacsgrp.c
	getacskeys.c
	getccdtab.c
	getcorner.c
	getlt.c
	getrefname.c
	getswitch.c
	gotfilename.c
	imghistory.c
	imgpedigree.c
	initmulti.c
	interp1d.c
	interpinfo.c
	key.c
	loadhead.c
	logprogress.c
	lowcase.c
	median.c
	mkname.c
	mkoutname.c
	mkspt.c
	mult1d.c
	multk1d.c
	multk2d.c
	omitstep.c
	parseamps.c
	parsedate.c
	prinfo.c
	reffiles.c
	rowpedigree.c
	sameint.c
	spline.c
	streqic.c
	sub1d.c
	tabhistory.c
	tabpedigree.c
	timestamp.c
	toelectrons.c
	trim.c
	ucalver.c
	ufilename.c
	unbinline.c
	whicherror.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME_INCLUDE_DIR}
	PUBLIC ${PROJECT_SOURCE_DIR}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)
