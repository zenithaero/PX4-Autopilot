add_subdirectory(sbgECom)
include_directories(sbgECom/src/ sbgECom/common/)

px4_add_module(
	MODULE modules__sbg
	MAIN sbg
	SRCS
		SBG.cpp
		SBG.hpp
	DEPENDS
		sbgECom
		px4_work_queue
	)
