px4_add_module(
	MODULE modules__zenith_control
	MAIN zenith_control
	SRCS
		ZenithControl.cpp
		ZenithControl.hpp
	DEPENDS
		px4_work_queue
	)
